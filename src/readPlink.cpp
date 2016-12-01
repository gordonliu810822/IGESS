
#include "readPlink.hpp"




// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


float normalCFD(float value)
{
    return 0.5 * erfc(-value * M_SQRT1_2);
}


/*get the positions of the identifiers in the column names*/
Col<int> getPositions(vector <string> fields, vector<string>& identifiers){

    uword size = identifiers.size();
    Col<int> pos(size);
    for(uword i = 0; i < size; i++){
        pos[i] = -1;
        for(uword j = 0; j < fields.size(); j++){
            if(fields[j].compare(identifiers[i]) == 0){
                pos[i] = (int)j ;
                break;
            }
        }

    }
    //    for(int i = (int)size - 1; i >= 0; i--){
    //        if( pos[i] == -1){
    //            identifiers.erase(identifiers.begin() + i);
    //        }
    //    }
    return pos;

}


GenoInfo::GenoInfo(string stringname) {
    string famfile = stringname;
    famfile += ".fam";
    string bimfile = stringname;
    bimfile += ".bim";

    this -> N =  getLineNum(famfile);
    this -> P =  getLineNum(bimfile);

    Chroms chroms(bimfile, P);
    fvec y = read_phenotypes(famfile, N);

    unsigned* X = new unsigned[ N * P];
    clock_t t1 = clock();
    readPlink(stringname,N, P, X);
    cout<<"Finish Reading Plink file, time elapsed:"<< (clock() - t1)*1.0/CLOCKS_PER_SEC << endl;
    cout <<"Sample Size =" <<  N << " SNP Number:" << P << endl;
    arma::Mat<unsigned>* Xdata = new arma::Mat<unsigned>(X, N, P, false,false);
    // arma::fmat corMat = calCorr(Xdata, index, avbIndex, bandwidth);
    this -> X = *Xdata;
    this -> X.replace(3, 0);
    this -> chroms = chroms;
    this -> y = y;
    delete[] X;

}




int snps_overlap(vector<SNP>& chrom_x_i, vector<SNP>& chrom_s_i, Col<uword>& xindex, Col<uword>& sindex){

    vector<SNP> common_snp_in_x;
    vector<SNP> common_snp_in_ss;
    //  cout << chrom_x_i[0] << "  " << chrom_x_i[1] << endl;
    sort(chrom_s_i.begin(), chrom_s_i.end());
    sort(chrom_x_i.begin(), chrom_x_i.end());

    //  cout << chrom_x_i[0] << "  " << chrom_x_i[1] << endl;
    set_intersection(chrom_x_i.begin(),chrom_x_i.end(),chrom_s_i.begin(),chrom_s_i.end(),back_inserter(common_snp_in_x));

    set_intersection(chrom_s_i.begin(),chrom_s_i.end(), chrom_x_i.begin(),chrom_x_i.end(),back_inserter(common_snp_in_ss));

    //  xindex.resize(common_snp_in_x.size());
    Mat<uword> indexPair(common_snp_in_x.size(), 2);
    for(int i = 0; i < common_snp_in_x.size(); i++){
        indexPair(i,0) = common_snp_in_x[i].idx;
        indexPair(i,1) = common_snp_in_ss[i].idx;
    }

    uvec sorted_index = sort_index(indexPair.col(0));
    indexPair = indexPair.rows(sorted_index);
    xindex = indexPair.col(0);
    sindex = indexPair.col(1);

    return (int)common_snp_in_x.size();

}

fvec read_phenotypes(string filename, int N){
    std::ifstream ifs(filename.c_str());

    std::string line;
    fvec y(N);
    string snpname;
    float pheno = 0;
    string tmp_str;
    int idx = 0;
    int phenopos = 5;
    while(std::getline(ifs, line)) // read one line from ifs
    {
        std::istringstream iss(line); // access line as a stream
        vector<string> fields;
        boost::split( fields, line, boost::is_any_of(" "));

    //     chromsome = (int)atoi(fields[0].c_str());
        pheno = (int)atoi(fields[phenopos].c_str());
//        for(int i=0; i < fields.size(); i++){
//            cout << fields[i] <<"  ";
//        }
//        cout << pheno << endl;
        y[idx] = pheno;
        
        idx++;

    }
   // cout << y << endl;
    ifs.close();
    return y;
}

Chroms read_snpnames(string filename, int P){
    std::ifstream ifs(filename.c_str());

    std::string line;
    Chroms chroms(P);//snps(P);
    int chromsome;
    string snpname;
    vector <string> fields;
    clock_t t1 = clock();
    int i = 0;
    while(std::getline(ifs, line)) // read one line from ifs
    {
        std::istringstream iss(line); // access line as a stream
        boost::split( fields, line, boost::is_any_of(" \t *"));
        chromsome = (int)atoi(fields[0].c_str());
        snpname = fields[1];
        int a1 = *fields[4].c_str();
        int a2 = *fields[5].c_str();
        chroms.chromsome[i] = chromsome;
        SNP snp(snpname, (int)i, from_x);
        chroms.snps[i] = snp;
        chroms.A1Array[i] = a1;
        chroms.A2Array[i] = a2;
        i++;
        //  chroms.push(snpname, chromsome, a1, a2);

    }
    ifs.close();
    return chroms;
}

int chroms_overlap(Chroms& chrom_x, Chroms& chrom_s,Col<uword>& xindex, Col<uword>& sindex){
    int total_overlap = 0;

    vector<SNP> geno_snps(chrom_x.snps, chrom_x.snps + chrom_x.P);
    vector<SNP> summary_snps(chrom_s.snps, chrom_s.snps + chrom_s.P);

    total_overlap = snps_overlap(geno_snps, summary_snps,  xindex, sindex);
    return total_overlap;

}

#define keep_cols 0
#define keep_rows 1

template <typename T>
void keep_indices(arma::Mat<T>& mat_obj, int type, uvec index)  {
    arma::Mat<T> mat_tmp;
    if(type == keep_cols){
       mat_tmp = mat_obj.cols(index);
    }else{
       mat_tmp = mat_obj.rows(index);
    }
    mat_obj.reset();
    mat_obj = mat_tmp;
}

void Summary::cal_overlap(GenoInfo& genoinfo)
{
    Col<uword> xindex;
    Col<uword> sindex;

    Col<int> direction;
    direction.reset();

    chroms_overlap(genoinfo.chroms, *this -> chroms, xindex, sindex);
    keep_indices(genoinfo.X, keep_cols, xindex);
    keep_indices((*this -> lpsummary), keep_rows, sindex);


    Chroms  chrom((int)xindex.size());
    chrom.chromsome = genoinfo.chroms.chromsome.elem(xindex);
    chrom.A1Array = genoinfo.chroms.A1Array.elem(xindex);
    chrom.A2Array = genoinfo.chroms.A2Array.elem(xindex);
    for(int i = 0; i < chrom.P; i++){
        chrom.snps[i] = genoinfo.chroms.snps[xindex[i]];
    }
    genoinfo.chroms.clear();
    genoinfo.chroms = chrom;
    genoinfo.xindex = xindex;
    genoinfo.N = (int)genoinfo.X.n_rows;
    genoinfo.P = (int)genoinfo.X.n_cols;

  //  return genoobj;

}




Summary::Summary(string summaryfile, vector<string> identifiers, string snpidentifier,
                 string chromidentifier){
    this -> chrom_no = 1;
    this -> type = type;
    cout << "loading summary data...." << endl;
    cout << "summaryfile =" << summaryfile << endl;
    this -> P = getLineNum(summaryfile) - 1;
    std::ifstream ifs(summaryfile.c_str());
    std::string line;
    std::getline(ifs, line);
    vector <string> fields;
    boost::split( fields, line, boost::is_any_of(" \t *"));
    
    identifiers.push_back(chromidentifier);
    identifiers.push_back(snpidentifier);

    Col<int> pos = getPositions(fields, identifiers);
    clock_t t1 = clock();
    uword snp_index = pos[pos.size()-1];
    // identifiers.erase(identifiers.end() - 1);
    int chr_index = pos[pos.size()-2];
    this -> chroms = new Chroms(P);
    if(chr_index == -1){
        this -> chroms -> chrom_no = 1;
    }
    
    pos.set_size(pos.size()-2);
    gwasidentifiers = identifiers;
    lpsummary = new Mat<double>(P, pos.size(), fill::randn);
    for(uword i = 0; i < P; i++)
    {
        std::getline(ifs, line);
        boost::split( fields, line, boost::is_any_of(" \t *"));
        int chrom_id = 1;
        if(chr_index != -1){
            chrom_id = (int)atoi(fields[chr_index].c_str());
        }

        this -> chroms -> chromsome[i] = chrom_id;
        SNP snp(fields[snp_index], (int)i, from_ss);

        this -> chroms -> snps[i] = snp;
        for(uword j = 0; j < pos.size(); j++){
            if(pos[j] == -1) continue;
            string value = fields[pos[j]];
            if(value != "NA"){
                float v = (float)atof(value.c_str());
                (*lpsummary)(i,j) = v;
            }

        }
    }


    if ( chr_index != -1){
        uvec indices = find_unique(this -> chroms -> chromsome);
        this -> chrom_no = (int)indices.size();
    }
    this -> convert();


    cout <<"Read summary time is " << (clock() - t1)/CLOCKS_PER_SEC << endl;
}




bool SNP::operator<(const SNP& obj) const{
    return this -> name < obj.name;
}
bool SNP::operator>(const SNP& obj) const{
    return this -> name > obj.name;
}

bool SNP::operator != (const SNP& obj) const{
    return this -> name.compare(obj.name) > 0;
}


bool SNP::operator == (const SNP& obj) const{
    return this -> name.compare(obj.name) == 0;
}



//GenoInfo ReadDataFromFile(string stringname) {
//  string famfile = stringname;
//  famfile += ".fam";
//  string bimfile = stringname;
//  bimfile += ".bim";
//
//  int N =  getLineNum(famfile);
//  int P =  getLineNum(bimfile);
//
//  Chroms chroms = read_snpnames(bimfile, P);
//  fvec y = read_phenotypes(famfile, N);
//
//  arma::vec index(22);
//  double sum2 = 0;
//  for(int i=1; i <= 22; i++ ){
//    index[i-1] = chroms.chromsomes[i-1].size(); //sum(*snps.chromsome == i);
//    sum2 += index[i-1];
//    cout <<"number of snps in chromsome "<<i <<" =" << index[i-1] << endl;
//  }
//  unsigned* X = new unsigned[ N * P];
//  clock_t t1 = clock();
//  readPlink(stringname,N, P, X);
//  cout<<"Finish Reading Plink file, time elapsed:"<< (clock() - t1)*1.0/CLOCKS_PER_SEC << endl;
//  cout <<"Sample Size =" <<  N << " SNP Number:" << P << endl;
//  arma::Mat<unsigned>* Xdata = new arma::Mat<unsigned>(X, N, P, false,false);
// // arma::fmat corMat = calCorr(Xdata, index, avbIndex, bandwidth);
//  GenoInfo obj;
//  obj.X = *Xdata;
//  obj.chroms = chroms;
//  obj.index = index;
//  obj.y = y;
//
//  return obj;
//}
