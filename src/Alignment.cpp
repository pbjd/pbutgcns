#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <cassert>
#include "Alignment.hpp"

using namespace dagcon;

Alignment::Alignment() : 
    start(0), 
    end(0), 
    frgid(""), 
    qstr(""), 
    tstr("") { }


void parseUnitig(std::istream& stream, Unitig* utg) {
    std::string line;
    std::getline(stream, line);
    std::stringstream row(line);
    std::string col;
    std::vector<std::string> fields;

    while(std::getline(row, col, ' ')) {
        if (col == "") continue;
        fields.push_back(col);
    }

    utg->id = fields[0];
    utg->seq = fields[1];
}

void parseLayout(std::istream& stream, Alignment* aln) {
    std::string line;
    std::getline(stream, line);
    std::stringstream row(line);
    std::string col;
    std::vector<std::string> fields;
    while(std::getline(row, col, ' ')) {
        if (col == "") continue;
        fields.push_back(col);
    }

    if (fields.size() == 0) return;

    //  tid, strand, tlen, tstart, tend, qstr, tstr
    aln->frgid = fields[0];

    std::istringstream ssStart(fields[1]);
    ssStart >> aln->start;

    std::istringstream ssEnd(fields[2]);
    ssEnd >> aln->end;

    aln->qstr = fields[3];
}

std::istream& operator>>(std::istream& instrm, Alignment& data) {
    parseLayout(instrm, &data);
    return instrm;
}

Alignment normalizeGaps(Alignment& aln) {
    size_t qlen = aln.qstr.length(), tlen = aln.tstr.length();
    assert(qlen == tlen);
    std::string qNorm, tNorm;

    // convert mismatches to indels
    for (size_t i=0; i < qlen; i++) {
        char qb = aln.qstr[i], tb = aln.tstr[i];
        if (qb != tb && qb != '-' && tb != '-') {
            qNorm += '-';
            qNorm += qb;
            tNorm += tb;
            tNorm += '-';
        } else {
            qNorm += qb;
            tNorm += tb;
        }
    }

    // update lengths
    qlen = qNorm.length();
    tlen = tNorm.length();

    // push gaps to the right, but not past the end
    for (size_t i=0; i < qlen-1; i++) {
        // pushing target gaps
        if (tNorm[i] == '-') {
            size_t j = i;
            while (true) {
                char c = tNorm[++j];
                if (c != '-' || j > qlen - 1) {
                    if (c == qNorm[i]) {
                        tNorm[i] = c;
                        tNorm[j] = '-';
                    }
                    break;
                }
            }
        }

        // pushing query gaps
        if (qNorm[i] == '-') {
            size_t j = i;
            while (true) {
                char c = qNorm[++j];
                if (c != '-' || j > tlen - 1) {
                    if (c == tNorm[i]) {
                        qNorm[i] = c;
                        qNorm[j] = '-';
                    }
                    break;
                }
            }
        }
    }

    // generate the final, normalized alignment strings
    Alignment finalNorm;
    finalNorm.frgid = aln.frgid;
    finalNorm.start = aln.start;
    for (size_t i=0; i < qlen; i++) {
        if (qNorm[i] != '-' || tNorm[i] != '-') {
            finalNorm.qstr += qNorm[i];
            finalNorm.tstr += tNorm[i];
        }
    }

    return finalNorm;
}
