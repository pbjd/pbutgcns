#include <cstdint>
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <fstream>
#include <log4cpp/Category.hh>
#include <log4cpp/Appender.hh>
#include <log4cpp/OstreamAppender.hh>
#include <log4cpp/Layout.hh>
#include <log4cpp/PatternLayout.hh>
#include <log4cpp/Priority.hh>
#include <boost/format.hpp>
#include <boost/call_traits.hpp>
#include <boost/progress.hpp>
#include <boost/bind.hpp>
#include <boost/program_options.hpp>
#include <boost/thread/thread.hpp>
#include "Alignment.hpp"
#include "AlnGraphBoost.hpp"
#include "SimpleAligner.hpp"
#include "BoundedBuffer.hpp"

namespace opts = boost::program_options;

typedef BoundedBuffer<dagcon::Alignment> AlnBuf;
typedef BoundedBuffer<dagcon::Alignment> CnsBuf;

void processUnitig(std::string fpath) {
    log4cpp::Category& logger = 
        log4cpp::Category::getInstance("processUnitig");

    std::ifstream ifs;
    ifs.open(fpath);

    Unitig utg;
    parseUnitig(ifs, &utg);
    logger.info("%s", utg.id.c_str());

    AlnGraphBoost ag(utg.seq);
    logger.info("Graph initialized. Length %d", utg.seq.size());

    logger.info("Aligning fragments");
    dagcon::Alignment aln;
    SimpleAligner align;
    int frgcount = 0;
    while (ifs >> aln) {
        aln.tstr = utg.seq.substr(aln.start, aln.end-aln.start);
        logger.debug("Aligning frg: %s, len: %d", 
            aln.frgid.c_str(), aln.qstr.size());
        align(aln);
        dagcon::Alignment norm = normalizeGaps(aln);
        ag.addAln(norm);
        ++frgcount;
    }
        
    logger.info("%d fragments aligned and added to graph", frgcount);

    logger.info("Generating consensus");
    ag.mergeNodes();
    std::string cns;
    ag.consensus(cns);
    boost::format fasta(">%s\n%s\n");
    fasta % utg.id; 
    fasta % cns;
    std::cout << fasta;
}

class Reader {
    Unitig utg_;
    AlnBuf* alnBuf_;
    std::istream* istrm_;
    int nAlnThreads_;
public:
    Reader(Unitig& utg, AlnBuf* b, std::istream* is, int n) : 
        utg_(utg),
        alnBuf_(b), 
        istrm_(is),
        nAlnThreads_(n)
    { }

    void operator()() {
        log4cpp::Category& logger = 
            log4cpp::Category::getInstance("Reader");

        dagcon::Alignment aln;
        int frgcount = 0;
        while (*istrm_ >> aln) {
            aln.tstr = utg_.seq.substr(aln.start, aln.end-aln.start);
            alnBuf_->push(aln);
            ++frgcount;
        }

        logger.info("%d fragments aligned and added to graph", frgcount);

        // write out sentinals, one per consensus thread
        dagcon::Alignment sentinel;
        sentinel.frgid = "_s_";
        for (int i=0; i < nAlnThreads_; i++)
            alnBuf_->push(sentinel);
    }
};

class GenAlign {
    AlnBuf* alnBuf_;
    CnsBuf* cnsBuf_;
public:
    GenAlign(AlnBuf* a, CnsBuf* b) :
        alnBuf_(a),
        cnsBuf_(b)
    { }

    void operator()() {
        log4cpp::Category& logger = 
            log4cpp::Category::getInstance("GenAlign");
        logger.info("Alignment thread started");
        SimpleAligner align;
        dagcon::Alignment aln;
        alnBuf_->pop(&aln);
        while (aln.frgid != "_s_") {
            logger.debug("Aligning frg: %s, len: %d", 
                aln.frgid.c_str(), aln.qstr.size());
            align(aln);
            cnsBuf_->push(normalizeGaps(aln));
            alnBuf_->pop(&aln);
        }
        // pass the sentinal
        cnsBuf_->push(aln);
    }
};

class Consensus {
    CnsBuf* cnsBuf_;
    Unitig draft_;
    int nAlnThreads_;
public:
    Consensus(CnsBuf* cp, const Unitig& utg, int n) : 
    cnsBuf_(cp),
    draft_(utg),
    nAlnThreads_(n)
    { }

    void operator()() {
        log4cpp::Category& logger = 
            log4cpp::Category::getInstance("Consensus");

        AlnGraphBoost ag(draft_.seq);
        logger.info("Graph initialized. Length %d", draft_.seq.size());

        dagcon::Alignment aln;
        int sentinelCount = 0;
        cnsBuf_->pop(&aln);
        while (true) {
            if (aln.frgid == "_s_" && ++sentinelCount == nAlnThreads_)
                break; 
            cnsBuf_->pop(&aln);
            ag.addAln(aln);
        }

        logger.info("Merging graph");
        ag.mergeNodes();

        logger.info("Generating consensus");
        std::string cns;
        ag.consensus(cns);

        boost::format fasta(">%s\n%s\n");
        fasta % draft_.id; 
        fasta % cns;
        std::cout << fasta;
    }

};

void setupLogger(log4cpp::Priority::Value priority) {
    // Setup the root logger to a file
    log4cpp::Appender *fapp = new log4cpp::OstreamAppender("default", &std::cerr);
    log4cpp::PatternLayout *layout = new log4cpp::PatternLayout();
    layout->setConversionPattern("%d %p [%c] %m%n");
    fapp->setLayout(layout); 
    log4cpp::Category& root = log4cpp::Category::getRoot();
    root.setPriority(priority);
    root.addAppender(fapp);
}

bool parseOpts(int ac, char* av[], opts::variables_map& vm) {
    opts::options_description 
        odesc("PacBio read-on-read error correction via consensus");
    odesc.add_options()
        ("help,h", "Display this help")
        ("verbose,v", "Increase logging verbosity")
        ("threads,j", opts::value<int>(), "Number of consensus threads to use")
        ("input", opts::value<std::string>(), "Unitig input file")
    ;

    opts::positional_options_description pdesc; 
    pdesc.add("input", 1);
    opts::store(opts::command_line_parser(ac, av)
                .options(odesc).positional(pdesc).run(), vm);

    opts::notify(vm);

    if (vm.count("help") || ! vm.count("input")) {
        std::cout << "Usage: " << av[0] << " [options] <input>"<< "\n\n";
        std::cout << odesc << "\n";
        return false;
    }


    return true;
}

int main(int argc, char* argv[]) {
    opts::variables_map vm;
    if (! parseOpts(argc, argv, vm)) exit(1);

    // http://log4cpp.sourceforge.net/api/classlog4cpp_1_1Priority.html
    // defaults to INFO.
    setupLogger(vm.count("verbose") ? 700 : 600);
    log4cpp::Category& logger = log4cpp::Category::getInstance("Main");
    std::string input = vm["input"].as<std::string>(); 
    if (vm.count("threads")) {
        int nthreads = vm["threads"].as<int>();
        logger.info("Multi-threaded. Input: %s, Threads: %d", 
            input.c_str(), nthreads);

        AlnBuf alnBuf(20);
        CnsBuf cnsBuf(20);

        std::ifstream ifs;
        ifs.open(input);
        Unitig utg;
        parseUnitig(ifs, &utg);
        logger.info("%s", utg.id.c_str());

        Consensus cns(&cnsBuf, utg, nthreads);
        boost::thread cnsThread(cns);

        std::vector<boost::thread> alnThreads;
        for (int i=0; i < nthreads; i++) {
            GenAlign ga(&alnBuf, &cnsBuf);
            alnThreads.push_back(boost::thread(ga));
        }

        Reader reader(utg, &alnBuf, &ifs, nthreads);
        boost::thread readerThread(reader);

        cnsThread.join();
        std::vector<boost::thread>::iterator it;
        for (it = alnThreads.begin(); it != alnThreads.end(); ++it)
            it->join();
    
        readerThread.join();
        
    } else {
        logger.info("Single-threaded. Input: %s", input.c_str());
        processUnitig(input);
    }

    return 0;
}
