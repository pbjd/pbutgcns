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
#include "Alignment.hpp"
#include "AlnGraphBoost.hpp"
#include "SimpleAligner.hpp"

namespace opts = boost::program_options;

void processUnitig(std::string fpath) {
    log4cpp::Category& logger = 
        log4cpp::Category::getInstance("processUnitig");

    std::ifstream ifs;
    ifs.open(fpath);

    Unitig utg;
    parseUnitig(ifs, &utg);
    logger.info("%s", utg.id.c_str());

    AlnGraphBoost ag(utg.seq);
    logger.info("Graph initialized. Adding fragments ...");

    dagcon::Alignment aln;
    SimpleAligner align;
    int frgcount = 0;
    while (ifs >> aln) {
        aln.tstr = utg.seq.substr(aln.start, aln.end-aln.start);
        align(aln);
        dagcon::Alignment norm = normalizeGaps(aln);
        ag.addAln(norm);
        ++frgcount;
    }
        
    logger.info("Loaded %d fragments", frgcount);
    ag.mergeNodes();

    std::vector<CnsResult> seqs;
    logger.info("Generating consensus", frgcount);
    ag.consensus(seqs, 0, 0);
    for (auto it = seqs.begin(); it != seqs.end(); ++it) {
        CnsResult result = *it;
        boost::format fasta(">%s/%d_%d\n%s\n");
        fasta % utg.id % result.range[0] % result.range[1];
        fasta % result.seq;
        std::cout << fasta;
    }
}

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
    std::string input = vm["input"].as<std::string>(); 
    processUnitig(input);

    return 0;
}
