#include <bits/stdc++.h>
using namespace std;

static const unordered_map<string, char> CODON_TABLE = {
    {"ATA",'I'},{"ATC",'I'},{"ATT",'I'},{"ATG",'M'},
    {"ACA",'T'},{"ACC",'T'},{"ACG",'T'},{"ACT",'T'},
    {"AAC",'N'},{"AAT",'N'},{"AAA",'K'},{"AAG",'K'},
    {"AGC",'S'},{"AGT",'S'},{"AGA",'R'},{"AGG",'R'},
    {"CTA",'L'},{"CTC",'L'},{"CTG",'L'},{"CTT",'L'},
    {"CCA",'P'},{"CCC",'P'},{"CCG",'P'},{"CCT",'P'},
    {"CAC",'H'},{"CAT",'H'},{"CAA",'Q'},{"CAG",'Q'},
    {"CGA",'R'},{"CGC",'R'},{"CGG",'R'},{"CGT",'R'},
    {"GTA",'V'},{"GTC",'V'},{"GTG",'V'},{"GTT",'V'},
    {"GCA",'A'},{"GCC",'A'},{"GCG",'A'},{"GCT",'A'},
    {"GAC",'D'},{"GAT",'D'},{"GAA",'E'},{"GAG",'E'},
    {"GGA",'G'},{"GGC",'G'},{"GGG",'G'},{"GGT",'G'},
    {"TCA",'S'},{"TCC",'S'},{"TCG",'S'},{"TCT",'S'},
    {"TTC",'F'},{"TTT",'F'},{"TTA",'L'},{"TTG",'L'},
    {"TAC",'Y'},{"TAT",'Y'},{"TAA",'*'},{"TAG",'*'},
    {"TGC",'C'},{"TGT",'C'},{"TGA",'*'},{"TGG",'W'}
};

string toupper_str(const string &s){
    string r=s; for(char &c:r) c = toupper((unsigned char)c); return r;
}
string revcomp(const string &s){
    string r; r.reserve(s.size());
    for(int i=(int)s.size()-1;i>=0;--i){
        char c = toupper((unsigned char)s[i]);
        char rc = 'N';
        if(c=='A') rc='T';
        else if(c=='T') rc='A';
        else if(c=='C') rc='G';
        else if(c=='G') rc='C';
        r.push_back(rc);
    }
    return r;
}
string translate(const string &dna){
    string aa; aa.reserve(dna.size()/3+1);
    for(size_t i=0;i+2<dna.size();i+=3){
        string cod = dna.substr(i,3);
        auto it = CODON_TABLE.find(cod);
        if(it!=CODON_TABLE.end()) aa.push_back(it->second);
        else aa.push_back('X');
    }
    return aa;
}

struct ORF { int start, end; bool forward; int frame; string nuc, aa; };

vector<ORF> find_orfs_forward(const string &s, int min_len){
    vector<ORF> res;
    int N = (int)s.size();
    for(int frame=0;frame<3;++frame){
        for(int i=frame;i+2<N;++i){
            if(s[i]=='A' && s[i+1]=='T' && s[i+2]=='G'){ // ATG
                int j=i+3;
                bool found=false;
                for(; j+2<N; j+=3){
                    string cod = s.substr(j,3);
                    if(cod=="TAA"||cod=="TAG"||cod=="TGA"){ found=true; break; }
                }
                if(found){
                    int start1 = i+1;
                    int end1 = j+3; // inclusive 1-based
                    int len = end1 - start1 + 1;
                    if(len>=min_len){
                        string nuc = s.substr(i, end1 - i);
                        string aa = translate(nuc);
                        res.push_back({start1, end1, true, frame, nuc, aa});
                    }
                    i = j; // continue after stop
                } else break;
            }
        }
    }
    return res;
}

vector<ORF> find_orfs_both(const string &orig, int min_len){
    string s = toupper_str(orig);
    string filtered;
    for(char c: s) if(c=='A'||c=='T'||c=='G'||c=='C'||c=='N') filtered.push_back(c);
    if(filtered.empty()) return {};
    auto fwd = find_orfs_forward(filtered, min_len);
    // reverse complement ORFs and map coords back to original
    string rc = revcomp(filtered);
    auto orfs_rc = find_orfs_forward(rc, min_len);
    vector<ORF> mapped_rc;
    int N = (int)filtered.size();
    for(auto &o: orfs_rc){
        int rc_start = o.start; // 1-based in rc
        int rc_end = o.end;     // 1-based in rc
        // map back: original positions = N - rc_end +1 ... N - rc_start +1
        int orig_start = N - rc_end + 1;
        int orig_end   = N - rc_start + 1;
        ORF mo = {orig_start, orig_end, false, o.frame, o.nuc, o.aa};
        mapped_rc.push_back(mo);
    }
    // combine and return
    vector<ORF> all;
    all.insert(all.end(), fwd.begin(), fwd.end());
    all.insert(all.end(), mapped_rc.begin(), mapped_rc.end());
    return all;
}

void print_orfs(const vector<ORF>& orfs){
    if(orfs.empty()){ cout << "No ORFs found with given minimum length.\n"; return;}
    cout << "Found " << orfs.size() << " ORFs:\n";
    for(size_t i=0;i<orfs.size();++i){
        const ORF &o = orfs[i];
        cout << "\nORF_" << (i+1) << (o.forward ? " (forward)":" (reverse)") << "\n";
        cout << "  Frame: " << o.frame << "\n";
        cout << "  Start (1-based): " << o.start << "\n";
        cout << "  End   (1-based): " << o.end << "\n";
        cout << "  Nuc len: " << (o.end - o.start + 1) << " nt\n";
        cout << "  AA len : " << o.aa.size() << " aa\n";
        cout << "  Nuc seq: " << o.nuc << "\n";
        cout << "  AA seq : " << o.aa << "\n";
    }
}

int main(int argc, char **argv){
    int min_len = 150;
    string seqarg;
    for(int i=1;i<argc;++i){
        string a = argv[i];
        if(a=="-m" && i+1<argc){ min_len = stoi(argv[++i]); }
        else if(a=="-h"||a=="--help"){ cout<<"Usage: "<<argv[0]<<" [-m min_nt] [sequence]\n"; return 0;}
        else seqarg = a;
    }

    string seq;
    if(!seqarg.empty()){
        seq = seqarg;
    } else {
        cout << "Enter DNA sequence (A/T/G/C/N) or paste FASTA/plain sequence, then press Enter:\n";
        if(!getline(cin, seq)) { cerr<<"No input received. Exiting.\n"; return 1; }
        // If more lines in stdin, read them too (useful when piping)
        string rest;
        while(getline(cin, rest)) { if(!rest.empty()) seq += rest; else break; }
    }

    // Clean and uppercase
    string s = toupper_str(seq);
    string filtered;
    for(char c: s) if(c=='A'||c=='T'||c=='G'||c=='C'||c=='N') filtered.push_back(c);
    if(filtered.empty()){ cerr << "No valid nucleotides found (A/T/G/C/N). Exiting.\n"; return 1; }
    cout << "Sequence length: " << filtered.size() << " nt\n";
    cout << "Minimum ORF length: " << min_len << " nt\n";

    auto all_orfs = find_orfs_both(filtered, min_len);

    // For better readability sort by start
    sort(all_orfs.begin(), all_orfs.end(), [](const ORF &a, const ORF &b){
        if(a.start != b.start) return a.start < b.start;
        return a.forward && !b.forward;
    });

    print_orfs(all_orfs);
    return 0;
}
