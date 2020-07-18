#include "fonction.h"

using namespace std;


/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function to check the number of threads
// Args: string afeter the -t option
void check_nb_thread(unsigned int nb_thread)
{
    #ifdef _OPENMP
        if (nb_thread > omp_get_max_threads() | nb_thread <= 0)
    {
        cout << "Incorrect thread number" << endl;
	exit(20);
    }
    #endif
}




/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function to check the number of threads
// Args: string after the -b option
void check_nb_boot(unsigned int nb_boot)
{
    if ( (nb_boot <= 0))
    {
        cout<< "Error: negative or null bootstrap number " << endl;
	exit(21);
    }
}


/****************************************************************************************************************************/
/****************************************************************************************************************************/
// Function to create liste of transcripts
// Args : path to the transcript file
int Isoforms_couting(std::string& TranscriptFilePathName,const std::string& TmpDir,int boot,const std::string& output_name)
{
    ifstream Transcriptfile(TranscriptFilePathName.c_str());

    if (!Transcriptfile)
    {
        cerr << " No transcript file: " << TranscriptFilePathName <<endl;
	supp_fichier(TmpDir,boot,output_name);
	exit(22);
    }
    else
    {
        int Isoforms_number=0;
        vector<std::string> GeneList;
        std::map<std::string,int> Map_GeneID_nbeIsoformes;
        string LineTranscript;
        while(getline(Transcriptfile, LineTranscript))
        {
            stringstream line(LineTranscript); // Turn the string into a stream.
            std::string SeqName;
            getline(line, SeqName, '\t'); // get sequence name (>ENSP00000375710|Homo_sapiens,>ENSP00000375706|Homo_sapiens,ect)
//             TranscriptsList.push_back(SeqName);
            
            std::string GeneID;
            getline(line, GeneID, '\t'); // get gene name (1,1,1,..) or (ENSG00000204673,ENSG00000204673,...)
            if (std::find(GeneList.begin(), GeneList.end(), GeneID) == GeneList.end()) // GeneID is not in GeneList
            {
                GeneList.push_back(GeneID);
                Map_GeneID_nbeIsoformes[GeneID]=1;
            }
            else
            {
                Map_GeneID_nbeIsoformes[GeneID]+=1;
            }
        }

//         for (std::map<std::string, int >::iterator ii = MapGeneTrans.begin(); ii != MapGeneTrans.end(); ii++)
        for (std::map<std::string,int >::iterator ii = Map_GeneID_nbeIsoformes.begin(); ii != Map_GeneID_nbeIsoformes.end(); ii++) // to parse the keys
        {
            int IsoNum_n=Map_GeneID_nbeIsoformes[ii->first];
            if (IsoNum_n!=1)
            {
                Isoforms_number+=IsoNum_n;
            }
        }
        return Isoforms_number;
    }
}
   

/****************************************************************************************************************************/
/****************************************************************************************************************************/		    
// Function  to check the input isoforms locus tag file 
// Args : path to the transcript file, path to the infile, temporary directory, number of bootstrap, the name of output
void checkTranscriptFile(std::string& TranscriptFilePathName,std::string& InfilePath, const std::string& TmpDir,int boot,const std::string& output_name)
{
    std::string InputTranscriptCMD="cut -f1 "+TranscriptFilePathName+" > "+TmpDir+"/Asupp/transcript_line_input.txt";
    system(InputTranscriptCMD.c_str());
    std::string file_line_number_input_str=TmpDir+"/Asupp/transcript_line_input.txt";
    ifstream file_line_number_input(file_line_number_input_str.c_str());
    if (!file_line_number_input)
    {
        cerr << " Bad format of transcript_file: " << TranscriptFilePathName <<endl;
        supp_fichier(TmpDir,boot,output_name);
        exit(23);
    }
    else
    {
        std::string line;
// 	     int nb_transcript=0;
        while(getline(file_line_number_input, line))
        {
            std::string grepCMDline="grep -F \""+line+"\" "+InfilePath+ " | wc -l";
            
            FILE* getRes = popen(grepCMDline.c_str(), "r");
            char buffer[128];
            std::string result = "";
            while (!feof(getRes)) 
            {
                if (fgets(buffer, 128, getRes) != NULL)
                    result += buffer;
            }
            
            if (result=="0\n")
            {
                cerr << "No transcript: \"" <<  line << "\" in the input dataset "<< endl << endl;
                supp_fichier(TmpDir,boot,output_name);
                exit(24);
            }
            pclose(getRes);
// 		nb_transcript+=1;
        }
// 	    return nb_transcript;
	}
}


/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function for fill the map CorrNumName
// Args: temporary directory
std::map<std::string,int> initMapCorNameNum(const std::string& TmpDir)
{
    string Corr_name_str = TmpDir+"/Asupp/Name_correspondance.txt";
    ifstream CorrName(Corr_name_str.c_str());

    if(CorrName)
    {
	std::string LineCorrName;
	std::map<std::string,int> map_name_num;
	while(getline(CorrName, LineCorrName))
	 {
	     stringstream msa(LineCorrName); // Turn the string into a stream.
	     std::string SeqName;
	     getline(msa, SeqName, '\t'); // get sequence name ">ENSP...."
	     
	     std::string SSeqNum;
	     getline(msa, SSeqNum, '\t'); // get the sequence number 1,2,5...
	     int SeqNum=0;
	     SeqNum=atoi(SSeqNum.c_str());
	     
	     map_name_num[SeqName]=SeqNum;
	 }
	 return map_name_num;
    }
    else
    {
	cerr << "can not open:" << Corr_name_str << endl;
    }
}


/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function for fill the map MapCorNumName
// Args: temporary directory
std::map<int, std::string > initMapCorNumName(const std::string& TmpDir)
{
   std::string Corr_name_str = TmpDir+"/Asupp/Name_correspondance.txt";
   ifstream CorrName(Corr_name_str.c_str());

    if(CorrName)
    {
	std::string LineCorrName;
	std::map<int, std::string >  map_numb_name;
	while(getline(CorrName, LineCorrName))
	 {
	     stringstream msa(LineCorrName); // Turn the string into a stream.
	     std::string SeqName;
	     getline(msa, SeqName, '\t'); // get sequence name ">ENSP...."
	     
	     std::string SSeqNum;
	     getline(msa, SSeqNum, '\t'); // get the sequence number 1,2,5...
	     int SeqNum=0;
	     SeqNum=atoi(SSeqNum.c_str());
	     
	     map_numb_name[SeqNum]=SeqName;
	 }
	 return map_numb_name;
    }
    else
    {
	cerr << "can not open:" << Corr_name_str << endl;
    }
}


/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function for fill the map TransNumSize
// Args: temporary directory
void initTransNumSize(const std::string& TmpDir,std::map<int,float>& map_name_num)
{
    string Corr_name_str = TmpDir+"/Asupp/Name_correspondance.txt";
    ifstream CorrName(Corr_name_str.c_str());

    if(CorrName)
    {
	std::string LineCorrName;
	while(getline(CorrName, LineCorrName))
	 {
	     stringstream msa(LineCorrName); // Turn the string into a stream.
	     std::string SeqName;
	     getline(msa, SeqName, '\t'); // get sequence name ">ENSP...."
	     
	     std::string SSeqNum;
	     getline(msa, SSeqNum, '\t'); // get the sequence number 1,2,5...
	     int SeqNum=0;
	     SeqNum=atoi(SSeqNum.c_str());
	     
	     std::string SSsize;
	     getline(msa, SSsize, '\t'); // get the sequence size 256, 276...
	     float SeqSize=0;
	     SeqSize=atof(SSsize.c_str());
	     
	     map_name_num[SeqNum]=SeqSize;
	 }
    }
    else
    {
	cerr << "can not open:" << Corr_name_str << endl;
    }
}



/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function to calcule the median of sequence sizes
// Args: temporary directory
float MedianComputation(std::vector<int>& TransNumList,std::map<int, float>& TransNumSize, int nb_seq,bool only_uniq_prot)
{
    std::vector<int> Vect_seq_sizes;
    if(only_uniq_prot)
	{
	    
	    for(int unsigned i=1;i< (nb_seq+1);i++)
		{
		    if (!(std::find(TransNumList.begin(), TransNumList.end(), i)!=TransNumList.end()) )// it is NOT an alternative transcript
		    {    
			Vect_seq_sizes.push_back(TransNumSize.find(i)->second);
		    }
		}
	}
    else
	{
	    for(int unsigned i=1;i< (nb_seq+1);i++)
		{
		   Vect_seq_sizes.push_back(TransNumSize.find(i)->second);
		}
	}
    
    sort(Vect_seq_sizes.begin(), Vect_seq_sizes.end());
    if(Vect_seq_sizes.size() % 2 == 0)
	{
	    return (Vect_seq_sizes[Vect_seq_sizes.size()/2 - 1] + Vect_seq_sizes[Vect_seq_sizes.size()/2]) / 2;
	}
    else
	{
	    return Vect_seq_sizes[Vect_seq_sizes.size()/2];
	}
}





/****************************************************************************************************************************/
/****************************************************************************************************************************/
// Function for fill the map initMapGeneTrans
// Args: temporary directory 
std::map<std::string, std::vector<int> > initMapGeneTrans(const std::string& TmpDir,std::string& transcript_File_Name,std::map<std::string,int>& MapCorNameNum)
{
    ifstream Transcriptfile(transcript_File_Name.c_str());
    if (Transcriptfile)
    {
	vector<int> TranscriptsList;
	string LineTranscript;
	std::map<std::string, std::vector<int> > map_gene_transcript;
	while(getline(Transcriptfile, LineTranscript))
	{
	    stringstream line(LineTranscript); // Turn the string into a stream.
	    std::string SeqName;
	    getline(line, SeqName, '\t'); // get sequence name (>ENSP00000375710|Homo_sapiens,>ENSP00000375706|Homo_sapiens,ect)
	    
	    int transNum;
	    std::map<std::string,int>::iterator it;

	    it = MapCorNameNum.find(SeqName);
	    if (it == MapCorNameNum.end())
	    {
		cerr << "Error in the MapCorNameNum key " << endl;
	    }
	    else
	    {
		transNum= MapCorNameNum.find(SeqName)->second;
		TranscriptsList.push_back(transNum);
		
		std::string Skey;
		getline(line, Skey, '\t'); // get gene name (1,1,1,..) or (ENSG00000204673,ENSG00000204673,...)
		
		map_gene_transcript[Skey].push_back(transNum);
	    }
	}
	Transcriptfile.close();
	return map_gene_transcript;
    }
    else
    {
	cerr << " Can not open transcript file: " << transcript_File_Name << endl;
    }
}
   
/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function for fill the map CorrNumName
// Args: temporary directory
std::map<int, std::vector<int> > initMapNumAltTrans(std::map<std::string, std::vector<int> >& MapGeneTrans)
{

    std::map<int, std::vector<int> >  MapNumAltTrans;
    vector<int> transcriptsKept;
    for (std::map<std::string, std::vector<int> >::iterator ii = MapGeneTrans.begin(); ii != MapGeneTrans.end(); ii++) // to parse the keys
    {
	vector<int> AltTrans=MapGeneTrans[ii->first];
        int number_of_isoforms;
        number_of_isoforms=0;
	for(int unsigned i=0;i<AltTrans.size();i++)
	{
	    MapNumAltTrans[AltTrans[i]]=AltTrans;
            number_of_isoforms++;
	}
	if (number_of_isoforms==1)
        { 
            MapNumAltTrans.erase(AltTrans[0]);
        }
    }
    return MapNumAltTrans;
}

/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function for fill the map CorrNumName
// Args: temporary directory
void PrintDistFile(std::map<int, std::string >& MapCorNumName,std::ofstream& DistFile,float& DistSum,int SeqNum,int SeqSize)
{
    std::string SeqName;
    std::map<int, std::string>::iterator itSeqName;
    itSeqName = MapCorNumName.find(SeqNum);
    if (itSeqName == MapCorNumName.end())
	{
	    cerr << "Error in the MapCorNumName key " << endl;
	}
    else
	{
	    SeqName=itSeqName->second;
	}
    DistFile << SeqName << "\t" << DistSum << "\t" << SeqSize << endl;
}

/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function for renaming sequences
// Args: input file name, temporary directory, working directory
int Seq_name_to_Number(std::string& infilePath, const std::string& TmpDir, const std::string& WorkDir,int boot,const std::string& output_name,bool DS)
{
    ifstream infile(infilePath.c_str());
    
        if (!infile)
        {
            cerr << "Cannot read input file: " << infilePath << endl;
	    supp_fichier(TmpDir,boot,output_name);
	    exit(25);
        }
        else
        {
	    string Msa_input_for_muscle = TmpDir+"/Asupp/Infile_number.fasta";
	    ofstream InfileNumber(Msa_input_for_muscle.c_str(), ios::out | ios::trunc);  
	    string Corr_name_str = TmpDir+"/Asupp/Name_correspondance.txt";
	    ofstream CorrName(Corr_name_str.c_str(), ios::out | ios::trunc);
	    
	    if(InfileNumber) 
	    {
		string line;
		int n=1;
		if(DS)
		{
		    int taille=0;
		    // check fasta format 
		    getline(infile, line);
		    if(line[0]=='>')
			{
			    InfileNumber << ">" << n  << endl;
			    CorrName << line << "\t" << n << "\t" ;
			    n++; 
			}
		    else
			{
			    cout << " Bad format of input file: must be in fasta format." << endl;
			    supp_fichier(TmpDir,boot,output_name);
			    exit(26);
			}
		    
		    // Then rename sequences
		    while(getline(infile, line))  
		    {
			if(line[0]=='>')
			{
				    CorrName << taille << endl;
				    CorrName << line << "\t" << n << "\t" ;
				    InfileNumber << ">" << n  << endl;
				taille=0;
				n++;
			}
			else
			{
			    taille+=line.size();
			    InfileNumber << line << endl;
			}
		    }
		    // for the last sequence
		    CorrName << taille << endl;
		    CorrName << line << "\t" << n << "\t" ;
		}
		
		else
		{
		    // check fasta format 
		    getline(infile, line);
		    if(line[0]=='>')
			{
			    InfileNumber << ">" << n  << endl;
			    CorrName << line << "\t" << n << endl;
			    n++;
			}
		    else
			{
			    cout << " Bad format of input file: must be in fasta format." << endl;
			    supp_fichier(TmpDir,boot,output_name);
			    exit(27);
			}
		    
		    // Then rename sequences
		    while(getline(infile, line))  
		    {
			if(line[0]=='>')
			{
			    InfileNumber << ">" << n  << endl;
			    CorrName << line << "\t" << n << endl;
			    n++;
			}
			else
			{
			    InfileNumber << line  << endl;
			}
		    }
		}
                InfileNumber.close(); 
		CorrName.close();
		return n-1;
	    }
	    else
	    {
		cerr << "Cannot create file: " << Msa_input_for_muscle << endl;
		supp_fichier(TmpDir,boot,output_name);
		exit(28);
	    }
	}
	infile.close();
}




/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function to reorder the alignment as in the input file
void BaseOrder(const std::string& TmpDir,const std::string& WorkDir,const std::string& output_name)
{
    string MsaFileName = TmpDir+"/Asupp/infile_align_desordered.fasta";
    ifstream MsaFile(MsaFileName.c_str());
      
      if (!MsaFile)
        {
            cerr << " No file: " << MsaFileName << endl;
        }
        else
        {
	    std::string Formatted_basal_MSA_path = TmpDir+"/Asupp/align_with_nb.fasta";
	    std::ofstream Formatted_basal_MSA(Formatted_basal_MSA_path.c_str(), ios::out | ios::trunc);  // creating the file "/Asupp/align_with_nb.fasta"

	    std::string Format_Msa_Path= TmpDir+"/format_msa/base_format";  
	    std::ofstream Format_Msa(Format_Msa_Path.c_str(), ios::out | ios::trunc);  //
	    
	    if (!Formatted_basal_MSA)
	    {
		cerr << " Cannot create file: " << Formatted_basal_MSA << endl;
	    }
        else
        {	    
	    if(!Format_Msa)
	    {
		cerr << " Cannot create file: " << Format_Msa << endl;
	    }
	    else
	    {
		std::string line;
		getline(MsaFile, line, '>'); // split the fasta file by sequences, first is an empty line
		std::map<int,std::string> my_map; // dictionnary of sequences name: sequence
		std::string nom_fichier;
		
		while(getline(MsaFile, line, '>'))  // split the fasta file by sequences
		{
		    stringstream msa(line); // Turn the string into a stream.
		    std::string Skey;
		    getline(msa, Skey); // get sequence name (63,56,25...)
		    int key=0;
		    key=atoi(Skey.c_str());
		    
		    std::vector<std::string> SeqTot; // sequence split by \n
		    std::string seqLine;
		    while(getline(msa, seqLine)) 
		    {
			SeqTot.push_back(seqLine); // sequence as a vector
		    }
		    
		    std::string SeqTotale; //sequence as a string
		    for(int i=0; i < SeqTot.size(); i++)
		    {
			SeqTotale +=SeqTot[i];
		    }
		    
		    my_map[key]=SeqTotale; // enter the sequence as a string in the dictionnary
		}

		for (std::map<int, std::string>::iterator i = my_map.begin(); i != my_map.end(); i++) // print the ordered msa in the file
		{
		    Formatted_basal_MSA << ">" << i->first << "\n";
		    Formatted_basal_MSA << i->second << "\n";
		    Format_Msa << i->second << "\n";
		}
		my_map.clear();
	    }
	}
	}
	MsaFile.close();
}


/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function for renaming sequences from number to initial name
// Args: temporary directory, working directory and name of the output
void Seq_num_to_seq_names(const std::string& TmpDir,const std::string& WorkDir,const std::string& output_name,int boot)
{
     string Corr_name_str = TmpDir+"/Asupp/Name_correspondance.txt";
     ifstream CorrName(Corr_name_str.c_str());  //déclaration du flux et ouverture du fichier
     
     string Align_number_str = TmpDir+"/Asupp/align_with_nb.fasta";
     ifstream AlignNumb(Align_number_str.c_str());  //déclaration du flux et ouverture du fichier
     
     if (!CorrName)
     {
	 cout << "Cannot read: " << Corr_name_str << endl;
	 supp_fichier(TmpDir,boot,output_name);
	 exit(29);
     }
     else
     {
	 string LineCorrName;
	 std::map<int,std::string> map_numb_name;
	 std::string nom_fichier;
	 while(getline(CorrName, LineCorrName))
	 {
	     stringstream msa(LineCorrName); // Turn the string into a stream.
	     std::string SeqName;
	     getline(msa, SeqName, '\t'); // get sequence name (63,56,25...)

	     std::string Skey;
	     getline(msa, Skey, '\t');
	     int key=0;
	     key=atoi(Skey.c_str());
	     
	     map_numb_name[key]=SeqName;
	 }

	if (!AlignNumb)
	{
	    cout << "Cannot read: " << Align_number_str << endl;
	    supp_fichier(TmpDir,boot,output_name);
	    exit(30);
	} 
	
	string AlignNameStr = WorkDir+output_name+".aln";
	ofstream AlignName(AlignNameStr.c_str(), ios::out | ios::trunc); 
	
	if (!AlignName)
	{
	    cout << "Cannot create file: " << AlignNameStr << endl;
	    supp_fichier(TmpDir,boot,output_name);
	    exit(31);
	}
	else
	{ 
	    string LineAlignNumb;
	    while(getline(AlignNumb, LineAlignNumb))
	    {
		if(LineAlignNumb[0]=='>')
		{
		    int k;
		    std::string num_str =std::string(LineAlignNumb.begin() + 1, LineAlignNumb.end());
		    k=atoi(num_str.c_str());
		    AlignName << map_numb_name.find(k)->second << endl; 
		}
		else
		{
		    AlignName << LineAlignNumb  << endl;
		}
	    }
	    AlignName.close();
	}
	map_numb_name.clear();
     }
     CorrName.close();
     AlignNumb.close();
}
 
  
/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Return the list of the transcript numbers
// Args: the Map containg the gene name and the transcripts number
std::vector<int> ListTranscriptNumber(std::map<int, std::vector<int> >& MapNumAltTrans)
{
     std::vector<int> ListTrans;
     std::map<int, std::vector<int> >::iterator it;
     for(it = MapNumAltTrans.begin();it != MapNumAltTrans.end(); ++it)
     {
	 int TransNum;
	 TransNum=it->first; // TransNum = 1 or 2 or 3 or.... n
	 ListTrans.push_back(TransNum);
     }
     return ListTrans;
 }
 
 
 
  
/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Reading the reference msa
// Args: file of the reference alignment, temporary directory
std::vector<std::string> create_align_tab(ifstream& file)
{
    std::vector<std::string> tab;
    std::string ligne;

    while(getline(file, ligne))
	{
	    tab.push_back(ligne);
	}
    return tab;
}



 
/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Creating array with number as residues of reference MSA
// Args: alignment's array, number of seq in input dataset, number of sites in the reference MSA
void create_tab_num_base(float** tab_num, std::vector<std::string> tab, unsigned nb_seq, unsigned nb_site)
{
    int aa(1);
     // For each residue we assign a the position number and 0 if it is a gap
    for(size_t i(0) ; i < nb_seq ; ++i)
    {
	for(size_t j(0) ; j < nb_site ; ++j)
	{
	    if(tab[i][j] == '-' )
	    {
		tab_num[i][j] = 0;
	    }
	    else
	    {
		tab_num[i][j] = aa;
		++aa;
	    }
	}
	aa = 1;
    }
}

/****************************************************************************************************************************/
/****************************************************************************************************************************/

// If option -auto, this function determines the number of gaps per site
// Args: alignment's array, number of seq in input dataset, number of sites in the reference MSA
bool stats_alignment(float** tab_num, unsigned nb_seq, unsigned nb_site)
{
    #ifdef _OPENMP
    float gap_j(0);
    float bad_site(0);
    std::vector<int> tailles_seq(nb_seq,0);
    float gap_percent=0.0;
    float percent_bad_site=0.0;
    
 	for(size_t j(0) ; j < nb_site ; ++j) // loop on the sites
	{
	    for(size_t i(0) ; i < nb_seq ; ++i) // loop for the sequences
	    {
		if(tab_num[i][j] == 0 )
		{
		    gap_j++;
		}
		else
		{
		    if(tab_num[i][j] > tailles_seq[i])
		    {
			tailles_seq[i]=tab_num[i][j];
		    }
		}
	    }
	    gap_percent=(gap_j/(1.0*nb_seq));
	    if( gap_percent > 0.8 )
	    {
		bad_site++;
	    }
	    gap_j=0.0;
	}

	percent_bad_site=bad_site/(1.0*nb_site);

	if( (percent_bad_site) > 0.35)
	{
	    cout << " Use of the -gap option, there is " << percent_bad_site*100 << " % of sites with more than 80% of gaps " << endl;
	    return true;
	}
	else
	{
	    return false;
	}
    #endif
}

/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function to calcul the mean length of sequences
// Args: vector of sequences sizes, number of sequences
float MeanComputation(std::vector<int> tailles_seq, unsigned nb_seq)
{
    double sum=0;
    float mean=0;
    
    for(size_t i(0) ; i < nb_seq ; ++i) // loop on the seq
    {
	sum+=tailles_seq[i];
    }
    mean=1.0*(sum/nb_seq);
    
    return(mean);
}

/*****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function to create the array containing the number of residue at each position of the reference alignment
// Args: array to fill, reference alignment,number of seq in input dataset, number of sites in the reference MSA
void create_tab_nb_non_gap(float* tab_nb_non_gap, float** tab_num_base, unsigned nb_seq, unsigned nb_site_base)
{
    float nb = 0.0;
    for (unsigned site = 0 ; site < nb_site_base ; ++site)
    {
        nb = 0;
        for (unsigned seq = 0 ; seq < nb_seq ; ++seq)
        {
            if (tab_num_base[seq][site] != 0.0)
            {
                ++nb;
            }
        }

        tab_nb_non_gap[site] = nb;
    }
}


/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Rename a file for fastdist
// Argument: the file, the temporary directory
void rename_for_fastdist(ifstream& base_MSA_file, const std::string& TmpDir,bool DS)
{
        if (!base_MSA_file)
        {
            cerr << "No file \"_base_msa.fasta\" " << endl;
        }
        else
        {
	    if(DS)
	    {
		string ResultFileFastdist = TmpDir+"/Asupp/align_with_nb_gap_X.fasta";
		ofstream pour_fastdist(ResultFileFastdist.c_str(), ios::out | ios::trunc);  
		if(pour_fastdist) 
		    {
			string line;
			int n=1;
			while(getline(base_MSA_file, line))
			{
			    if(line[0]=='>')
			    {
				pour_fastdist << ">" << n  << endl;
				n++;
			    }
			    else
			    {
				pour_fastdist << StrReplace(line,"-","X") << endl;
			    }
			}
			pour_fastdist.close();  
		    }
		else
		{
		    cerr << "Error opening file: " << TmpDir <<"/Asupp/align_with_nb_gap_X.fasta base_MSA_file"<< endl;
		}
		base_MSA_file.close();
	    }
	    else
	    {
		string ResultFileFastdist = TmpDir+"/Asupp/align_with_nb.fasta";
		ofstream pour_fastdist(ResultFileFastdist.c_str(), ios::out | ios::trunc);  
		if(pour_fastdist) 
		{
		    string line;
		    int n=1;
		    while(getline(base_MSA_file, line))
		    {
			if(line[0]=='>')
			{
			    pour_fastdist << ">" << n  << endl;
			    n++;
			}
			else
			{
			    pour_fastdist << line << endl;
			}
		    }
                pour_fastdist.close();  
		}
		else
		{
		    cerr << "Error opening file: " << TmpDir <<"/Asupp/align_with_nb_gap_X.fasta base_MSA_file"<< endl;
		}
		base_MSA_file.close();
	    }
    }
}





/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function to replace string in string
string StrReplace(std::string& str, const std::string& oldStr, const std::string& newStr)
{
  size_t pos = 0;
  while((pos = str.find(oldStr, pos)) != std::string::npos)
  {
     str.replace(pos, oldStr.length(), newStr);
     pos += newStr.length();
  }
  return str;
}


/****************************************************************************************************************************/
/****************************************************************************************************************************/

// If NAN in the distance matrices
// Args: the matrice file to modify, the temporary directory, the number of seq in input dataset, the bootstrap number in process
// void replace_NAN(std::ifstream& matrices_file, const std::string& TmpDir,int nb_seq,int boot,bool DS)
// {
//         if (!matrices_file)
//         {
//             cerr << "No matrice file. Error in Fastdist run. "  << endl;
//         }
//         else
//         {
// 	    string line;
// 	    int numSeq=0;
// 	    int nb_tree=0;
// 	    while(numSeq < (nb_seq+1)) // the first matrice is the one of the reference alignment 
// 		{
// 		    if(DS)
// 		    {
// 			std::string ResultFile = TmpDir + "/Asupp/matrice_DistancesScores";
// 			std::ofstream without_NAN(ResultFile.c_str(), ios::out | ios::trunc); 
// 			if(without_NAN) 
// 			{
// 			    while(numSeq < (nb_seq+1))
// 			    {
// 				getline(matrices_file, line);
// 				without_NAN << StrReplace(line,"nan","6") << endl;
// 				numSeq++;
// 			    }
// 			}
// 			else
// 			{
// 			    cerr << "Error opening matrice_sans_nan" << endl;
// 			}
// 		    }
// 		    else
// 		    {
// 			  getline(matrices_file, line);
// 			  numSeq++;
// 		    }
// 		}
// 	    numSeq=0;
//             // instructions
// 	    for(int nb_tree(0); nb_tree < boot ; nb_tree++)
// 		{
// 		    std::stringstream ResultFile_str;
// 		    ResultFile_str << TmpDir << "/Asupp/matrice_" << nb_tree << "_sans_nan";
// 		    string ResultFile = ResultFile_str.str();
// 		    ofstream without_NAN(ResultFile.c_str(), ios::out | ios::trunc); 
// 		    if(without_NAN) 
// 			{
// 			    while(numSeq < (nb_seq+1))
// 				{
// 				    getline(matrices_file, line);
// 				    without_NAN << StrReplace(line,"nan","6") << endl;
// 				    numSeq++;
// 				}
// 			without_NAN.close(); 
// 			}
// 		    else
// 			{
// 			    cerr << "Error opening matrice_nan_" << nb_tree << endl;
// 			}
// 		    numSeq=0;
// 		}
// 	}
// 	matrices_file.close();
// }




/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function to check if the input tree is ok (not branch lenght in scientific format, etc)
// Args: the tree to check
std::string checkTree(std::string& tree)
{
    while( tree.find("e-") != -1)
    {
	int pos=tree.find("e-");
	int Colon=pos;
	while ( tree.at(Colon)!= ':')
	{
	    Colon -= 1 ;
	}
	
	int bracket=pos+1;
	while ( (tree.at(bracket)!=',') && (tree.at(bracket)!=')') )
	{
	    bracket+=1;
	}
	int len;
	len= bracket-Colon;
	tree.replace(Colon+1,len-1,"0");
    }
    
    while (tree.find("-") != -1)
    {
	int posNegBeg=tree.find("-");
	int posNegEnd=posNegBeg+1;
	while ( (tree.at(posNegEnd)!=(',')) && (tree.at(posNegEnd)!=(')')) )
	{
	    posNegEnd+=1;
	}
	tree.replace(posNegBeg,posNegEnd-posNegBeg,"0.00000");
    }
   
    return tree;
}



/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function to find to sister leaf in a tree
// Args: the tree, the output tree
std::string FindSisterLeafs(std::string& tree,std::ofstream& RootedForMafft)
{
    int bracket1=tree.find('(');
    int bracket2=bracket1+1;
    int End=tree.size();
    while (bracket2!=End)
    {
	if (tree.at(bracket2)=='(')
	{
	    bracket1=bracket2;
	    bracket2+=1;
	}
	else
	    if (tree.at(bracket2)==')')
	    {
		break;
	    }
	    else
	    {
		bracket2+=1;
	    }
    }

    int PosNode1=bracket1+1;
    while (tree.at(PosNode1)!=':')
    {
	PosNode1+=1;
    }
    
    int PosLenght1=PosNode1+1;
      while (tree.at(PosLenght1)!=',')
    {
	PosLenght1+=1;
    }  
    
    int PosNode2=PosLenght1+1;
    while (tree.at(PosNode2)!=':')
    {
	PosNode2+=1;
    }
    
    int PosLenght2=PosNode2+1;
    while (tree.at(PosLenght2)!=')') 
    {
	PosLenght2+=1;
    }  
    
    int Node1=atoi(tree.substr(bracket1+1,PosNode1-bracket1-1).c_str());
    int Node2=atoi(tree.substr(PosLenght1+1,PosNode2-PosLenght1).c_str());
    string Lenght1=tree.substr(PosNode1+1,PosLenght1-PosNode1-1);
    string Lenght2=tree.substr(PosNode2+1,PosLenght2-PosNode2-1);   
   
   if (Node1 < Node2 )
   { 
       RootedForMafft << setw(5) << Node1 << setw(6) << Node2 << setw(11) << Lenght1 << setw(11) << Lenght2 << endl;
       tree.replace(bracket1,bracket2-bracket1+2,tree.substr(bracket1+1,PosNode1-bracket1));
   }
   else
   {
      RootedForMafft << setw(5) << Node2 << setw(6) << Node1 << setw(11) << Lenght2 << setw(11) << Lenght1 << endl;
      tree.replace(bracket1,bracket2-bracket1+2,tree.substr(PosLenght1+1,PosNode2-PosLenght1));
   }
   return tree;
}




/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function to formate the rooted tree in the input guide tree mafft format 
// Args:  the temporary directory, the number of tree in process
void MafftTreeFormating(const std::string& TmpDir,int nb_tree)
{ 
    std::stringstream NamesString;
    NamesString << TmpDir << "/Asupp/rooted_tree_" << nb_tree;
    std::string NameRootTree = NamesString.str();
    std::ifstream RootedTree(NameRootTree.c_str());

    if (!RootedTree)
	{
		cerr << "Cannot open the file: " << NameRootTree << endl;
	}
	else
	{
	    std::string Tree;
	    getline(RootedTree, Tree);
	   
	    std::string TreeOk;
	    TreeOk=checkTree(Tree);
	    
	    std::stringstream RootedForMafftNamesString;
	    RootedForMafftNamesString << TmpDir << "/Asupp/rooted_tree_mafft_" << nb_tree;
	    std::string RootedForMafftName=RootedForMafftNamesString.str();
	    std::ofstream RootedForMafft(RootedForMafftName.c_str(), ios::out | ios::trunc); 
	    
	    while (TreeOk.find('(')!= -1 )
	    {
		TreeOk=FindSisterLeafs(Tree,RootedForMafft);
	    }
	    RootedForMafft.close();
	    
	    remove(NameRootTree.c_str());
	}
	RootedTree.close();
}



/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function to re order the perturbed msa
// Args:  the temporary directory, the number of tree in process
void MsaOrder(const std::string& TmpDir,int nb_tree)
{
    stringstream TreeSStr;
    TreeSStr << nb_tree;
    string Tree_Str = TreeSStr.str();
    string MsaFileName = TmpDir+"/perturbed_msa/msa_"+Tree_Str+".fasta";
    ifstream MsaFile(MsaFileName.c_str());
      
      if (!MsaFile)
        {
            cerr << " No file: " << MsaFileName << endl;
        }
        else
        {
	    std::string Formatted_MSA_path = TmpDir+"/format_msa/msa_"+Tree_Str;
	    std::ofstream Formatted_MSA(Formatted_MSA_path.c_str(), ios::out | ios::trunc);  // openning of file msa_X

	    if(Formatted_MSA)
	    {
		std::string line;
		getline(MsaFile, line, '>'); // split the fasta file by sequences, first is an empty line
		std::map<int,std::string> my_map; // dictionnary of sequences name: sequence
		std::string nom_fichier;
		
		while(getline(MsaFile, line, '>'))  // split the fasta file by sequences
		{
		    stringstream msa(line); // Turn the string into a stream.
		    std::string Skey;
		    getline(msa, Skey); // get sequence name (63,56,25...)
		    int key=0;
		    key=atoi(Skey.c_str());
		    
		    std::vector<std::string> SeqTot; // sequence split by \n
		    std::string seqLine;
		    while(getline(msa, seqLine)) 
		    {
			SeqTot.push_back(seqLine); // sequence as a vector
		    }
		    
		    std::string SeqTotale; //sequence as a string
		    for(int i=0; i < SeqTot.size(); i++)
		    {
			SeqTotale +=SeqTot[i];
		    }
		    
		    my_map[key]=SeqTotale; // enter the sequence as a string in the dictionnary
		}

		for (std::map<int, std::string>::iterator i = my_map.begin(); i != my_map.end(); i++) // print the ordered msa in the file
		{
		    Formatted_MSA << i->second << "\n";
		}
		my_map.clear();
	    }
	}
	MsaFile.close();
}





/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function creating the table of line numbers for analysis
// Args: transcript file name, array to fill, number of sequences to analyze
std::vector<short>  create_species_line(std::string& TranscriptFilePathName,std::string& InfilePath, const std::string& TmpDir,int boot,const std::string& output_name)
{
    ifstream Transcriptfile(TranscriptFilePathName.c_str());

    if (!Transcriptfile)
    {
        cerr << " No transcript file: " << TranscriptFilePathName <<endl;
	supp_fichier(TmpDir,boot,output_name);
	exit(32);
    }
    else
    {

	string SpeciesCMD = "grep -F '>' "+InfilePath+" | grep -F \"$( (cut -f 1 "+TranscriptFilePathName+") )\" -n | awk -F':' '{print $1}' > "+TmpDir+"/Asupp/transcript_line.txt";
	system(SpeciesCMD.c_str());
	string file_line_number_str=TmpDir+"/Asupp/transcript_line.txt";
	ifstream file_line_number(file_line_number_str.c_str());
	if (!file_line_number)
	{
	    cerr << " Unable to find sequences of transcript file " << file_line_number_str <<endl;
	    supp_fichier(TmpDir,boot,output_name);
	    exit(33);
	}
	else
	{
	    std::vector<short> tab;
	    std::string ligne;
	    while(getline(file_line_number, ligne))
	    {
		short test=(short) atoi(ligne.c_str());
		tab.push_back(test-1); // because numerotation of grep begin at 1
	    }
	    file_line_number.close();
	return tab;
	}
    }
    Transcriptfile.close();
}





/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Array of perturbed alignment
// Argument: file of the perturbed msa, the array to fill
void create_align_tab_pert(std::ifstream& file, std::string* tab)
{
    std::string ligne;
    unsigned int i(0);
    
     // each line is stocked in the array
    while(getline(file, ligne))
	{
	    tab[i] = ligne;
	    ++i;
	}
}




/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function to create the array of perturbed msa with numbers instead of residues
// Args: array to fill,number of seq in input dataset, number of sites in the reference MSA
void create_tab_num_pert(short** tab_num, std::string* tab, unsigned nb_seq, unsigned nb_site)
{
    int nb_aa(1);
     // For each residue we assign a the position number and 0 if it is a gap
    for(size_t i(0) ; i < nb_seq ; ++i)
	{
	    for(size_t j(0) ; j < nb_site ; ++j)
	    {
		if(tab[i][j] == '-' )
		    {
			tab_num[i][j] = 0;
		    }
		else
		    {
			tab_num[i][j] = nb_aa;
			++nb_aa;
		    }
	    }
    nb_aa = 1;
	}
}




/****************************************************************************************************************************/
/****************************************************************************************************************************/

// The main function, computing the residue scores for one perturbed alignment
// Args: array of the reference alignment,  array of the perturbed alignment, number of seq in input dataset, number of sites in the reference MSA, array to fill, species for which we want a score, number of residue at each position of the reference msa, GapPenalty option statut
void create_residu_score_for1msa(float** tab_num_base,short** tab_num_pert, unsigned nb_seq, unsigned nb_site, float** tab_residu_score_for1msa, std::vector<short> seq_lines,float* tab_nb_non_gap,bool GapPen)
{

    vector<short>::iterator seq;
    short nb_ligne(0);

     // iteration one each sequence for which we want a score
    for(seq = seq_lines.begin(); seq!=seq_lines.end(); ++seq) // Pour chaque séquence de l'alignement de référence
    {
        for (unsigned site(0) ; site < nb_site ; ++site) // Pour chaque site en question de l'alignement de référence
        {

            // Si le site dans l'alignement de base est un gap on donne un score de 5000 dans la nouvelle table
            if (tab_num_base[*seq][site] == 0.0)
            {
                tab_residu_score_for1msa[nb_ligne][site] = 5000 ;
            }
            // Sinon on calcule son score
            else
            {
                // On cherche dans l'alignement perturbé le même résidu de l'alignement de référence que l'on est en train de traiter et on stocke son site dans la variable nb
                int nb(0);
                while (tab_num_base[*seq][site] != tab_num_pert[*seq][nb])
                {
                    ++nb;
                }


                float somme (0.0), nb_paire_comptee(0.0);
                for (unsigned paire(0) ; paire < nb_seq ; ++paire) // Pour chaque residu paire
                {

                    // Si le residu paire est égal à 0 (donc un gap) ou si il est différent que le résidu en question
//                     bool compare1((tab_num_base[paire][site] != 0) & (paire != *seq));
                    if(paire != *seq)
                    {
			 if (tab_num_base[paire][site] != 0 ) // si dans alignement de base ce n'est pas à gap à la sequence paire
			 {
			    // On compare la paire de residu de l'alignement de base avec celle de l'alignement perturbé
			    bool compare2((tab_num_base[*seq][site] == tab_num_pert[*seq][nb]) & (tab_num_base[paire][site] == tab_num_pert[paire][nb]));

			    // Si les deux résidus sont les mêmes on assigne un score de 1 à la paire
			    if (compare2)
			    {
				if(GapPen)
				{
				    if(nb_seq != tab_nb_non_gap[site])
				    {
					somme = somme + (1.0 / ( (nb_seq) - (tab_nb_non_gap[site]) ) ); // on divise par le nombre de gaps à cet endroit
				    }
				    else // pas de gaps a cette position
				    {
					somme = somme + 1.0 ;
				    }
				}
				else
				{
				    somme = somme + 1.0;
				}
			    }
			 }
                     }
		 }

		 if(tab_nb_non_gap[site] != 1) // Si le site de la sequence n'est pas que face à des gap
		    {
			tab_residu_score_for1msa[nb_ligne][site] = somme / (tab_nb_non_gap[site]-1); // le score c'est le nombre de paire egale divisé par le nombre de residu à cette position dans le base_align
		    }
		    else 
			{
			    if (GapPen)
			    {
				tab_residu_score_for1msa[nb_ligne][site] = -1.0;
			    }
			    else
			    {
				tab_residu_score_for1msa[nb_ligne][site] = 1.0;
			    }
			}
		}
        }
    nb_ligne++;
    }
}





/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function to compute score for one bootstrap
// Args: array with score by seqeunces, array with residue scores, number of sequences to analyze, number of sites in the reference MSA, sequences to analyze, array of number of residue at each position in the reference msa, ShortPenalty option
void create_tab_score_seq(float* tab_score_par_seq, float** tab_score_res, unsigned nb_species, unsigned nb_seq, unsigned nb_site, std::vector<short> species_lines, float* tab_nb_non_gap,bool ShortPen)
{
    for(unsigned species(0) ; species < nb_species ; ++species) // for each sequence we compute the average residue scores
    {
        // Initialisation variables
        double somme(0.0), moyenne(0.0), nb(0.0), score_taille_aln(0.0), score_taille_seq(0.0);
	for(size_t site(0) ; site < nb_site; ++site)
	    {
		if (tab_score_res[species][site] != 5000) // it is not a gap 
		    {
			somme = somme + tab_score_res[species][site];
			++nb;
		    }
	    }

	if (ShortPen)
	    {
		score_taille_aln = ( somme / nb_site); // to penalyze short sequences
		tab_score_par_seq[species]=score_taille_aln;
	    }
	else
	    {
		moyenne = somme / nb;
		tab_score_par_seq[species] = moyenne;
	    }
    }
}



/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function to create the .scores file
void create_seq_score_file_with_lines(float* tab_score_per_res, std::string& InfilePath, vector<short> species_lines, unsigned nb_species,std::string out_name, const std::string& WorkDir,const std::string& TmpDir)
{

    string cmd ("sed '/^>/!d' "+InfilePath+" >"+TmpDir+"/Asupp/align_name");
    system(cmd.c_str());

    vector<short>::iterator seq;

    for(seq = species_lines.begin(); seq!=species_lines.end(); ++seq)
    {

        stringstream s;
        s << *seq+1;
        string str = s.str();
        string cmd2 = "sed -n -e '"+str+"p' "+TmpDir+"/Asupp/align_name >>"+TmpDir+"/Asupp/lines_names";
        system(cmd2.c_str());

    }

    string lines_names_str=TmpDir+"/Asupp/lines_names";
    ifstream names_file(lines_names_str.c_str());
    
    if (!names_file)
    {
	cerr << "No file: " << lines_names_str << endl;
    }
    
    else
    {
	string nom_fichier;
        nom_fichier = WorkDir+out_name+".scores";
        std::ofstream fichier_score_final(nom_fichier.c_str(), ios::out | ios::trunc); 

        if(fichier_score_final)
        {
            fichier_score_final << "SEQUENCE_NAMES" << "\t" << "SEQUENCE SCORE" << endl;
            for(unsigned int i(0) ; i < nb_species ; ++i)
            {
                string name;
                getline(names_file,name);
                fichier_score_final << name << "\t" ;
                fichier_score_final << tab_score_per_res[i] << endl;
            }
            fichier_score_final.close();
        }
        else
        {
            cerr << "Cannot open the file: " << nom_fichier << endl;
        }

    }
    names_file.close();
}




/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function to create the filtered output file
// Args: name of the output, sequences to analyze, path to the .scores output file, path to the input file, working directory
void filtered_file(const std::string& output_name,std::string& TranscriptFilePathName,std::string& ScoreFilePath, std::string& infilePath,const std::string& WorkDirectory,const std::string& TmpDir,int boot)
{
    ifstream Transcriptfile(TranscriptFilePathName.c_str());
    ifstream Infile(infilePath.c_str());
    ifstream ScoreFile(ScoreFilePath.c_str());
    
    if(!Transcriptfile)
    {
	cerr << " No transcript file: " << TranscriptFilePathName <<endl;
	supp_fichier(TmpDir,boot,output_name);
	exit(34);
    }
    else
	if (!Infile)
	{
	    cerr << " Error in infile: " << infilePath <<endl;
	    supp_fichier(TmpDir,boot,output_name);
	    exit(35);
	}
	else
	    if(!ScoreFile)
	    {
		cerr << " Error in infile: " << ScoreFilePath <<endl;
		supp_fichier(TmpDir,boot,output_name);
		exit(36);
	    }
	    else
	    {
		string FilteredFileName = WorkDirectory+output_name+"_filtered.fasta";
		ofstream FilteredFile(FilteredFileName.c_str(), ios::out | ios::trunc); 
		if(!FilteredFile)
		{
		    cerr << "cant not create: " << FilteredFileName << endl;
		}
		else
		{
		    vector<std::string> TranscriptsList;
		    string LineTrasncript;
		    std::map<std::string, vector<std::string> > map_gene_transcript;
		    while(getline(Transcriptfile, LineTrasncript))
		    {
			stringstream line(LineTrasncript); // Turn the string into a stream.
			std::string SeqName;
			getline(line, SeqName, '\t'); // get sequence name (>ENSP00000375710|Homo_sapiens,>ENSP00000375706|Homo_sapiens,ect)
			TranscriptsList.push_back(SeqName);
			
			std::string Skey;
			getline(line, Skey, '\t'); // get gene name (1,1,1,..) or (ENSG00000204673,ENSG00000204673,...)
			
			map_gene_transcript[Skey].push_back(SeqName);
		    }
		    
		    // creation of sequence name/score
		    string LineNameScore;
		    std::map<std::string,float> map_seq_score;
		    getline(ScoreFile, LineNameScore);
		    while(getline(ScoreFile, LineNameScore))
		    {
			stringstream lineBis(LineNameScore); // Turn the string into a stream.
			std::string SeqName;
			getline(lineBis, SeqName, '\t'); // get sequence name (>ENSP00000375710|Homo_sapiens,>ENSP00000375706|Homo_sapiens,ect)
			std::string score;
			getline(lineBis, score, '\t'); // get sequence score (0.98,0.95,...)
			
			map_seq_score[SeqName]=atof(score.c_str());
		    }
		    
		    vector<std::string> transcriptsKept;
		    for (std::map<std::string, vector<std::string> >::iterator ii = map_gene_transcript.begin(); ii != map_gene_transcript.end(); ii++) 
		    {
			float score_max=0.0;
			string best_seq;
			for( vector<string>::const_iterator eptr=ii->second.begin(); eptr!=ii->second.end(); eptr++)
			{
			    float score_seq=map_seq_score.find(*eptr)->second;
			    if(score_max==0.0)
			    {
				best_seq=*eptr;
				score_max=score_seq;
			    }
			    else
				if (score_seq>score_max)
				{
				    best_seq=*eptr;
				    score_max=score_seq;
				}
			}
			transcriptsKept.push_back(best_seq);
		    }
		    
		    string FastSeq;
		    getline(Infile,FastSeq,'>');
		    while(getline(Infile,FastSeq,'>'))
		    {
			stringstream FastSeqSS(FastSeq);
			string SeqName;
			getline(FastSeqSS,SeqName);
			SeqName=">"+SeqName;
			if( (std::find(TranscriptsList.begin(),TranscriptsList.end(),SeqName)==TranscriptsList.end()) || (std::find(transcriptsKept.begin(),transcriptsKept.end(),SeqName)!=transcriptsKept.end()) )
			    {
				FilteredFile << ">" << FastSeq << endl;
			    }
		    }
		}
		Transcriptfile.close();
		Infile.close();
		ScoreFile.close();
	    }
}




/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function to remove temporary files
// Args: temporary directory path, bootstrap number and output name
void supp_fichier(const std::string& TmpDir,int boot,const std::string& output_name)
{
    std::string Asupp=TmpDir+"/Asupp/";
    std::string perturbed_msa=TmpDir+"/perturbed_msa/";
    std::string format_msa=TmpDir+"/format_msa/";

    std::string InfileS=Asupp+"Infile_number.fasta";
    std::ifstream InfileF(InfileS.c_str());    
    if(InfileF)
    {
	remove(InfileS.c_str());
    }
    
    std::string Name=Asupp+"Name_correspondance.txt";
    std::ifstream NameF(Name.c_str());    
    if(NameF)
    {
	remove(Name.c_str());
    }
   
    std::string align=Asupp+"align_name";
    std::ifstream alignF(align.c_str());    
    if(alignF)
    {
	remove(align.c_str());
    }
    
    std::string alignWithNumber=Asupp+"align_with_nb.fasta";
    std::ifstream alignWithNumberF(alignWithNumber.c_str());    
    if(alignWithNumberF)
    {
	remove(alignWithNumber.c_str());
    }
    
    std::string lines_names=Asupp+"lines_names";  
    std::ifstream lines_namesF(lines_names.c_str());    
    if(lines_namesF)
    {
	remove(lines_names.c_str());
    }
      
    std::string matrices=Asupp+"tree.dst";
    std::ifstream matricesF(matrices.c_str());    
    if(matricesF)
    {
	remove(matrices.c_str());
    }
    
    std::string trees=Asupp+"tree.phb";
    std::ifstream treesF(trees.c_str());    
    if(treesF)
    {
	remove(trees.c_str());
    }
    std::string infile_align=Asupp+"infile_align_desordered.fasta";
    std::ifstream infile_alignF(infile_align.c_str());    
    if(infile_alignF)
    {
	remove(infile_align.c_str());
    }
    
    std::string transcript_line=Asupp+"transcript_line.txt";
    std::ifstream transcript_lineF(transcript_line.c_str());    
    if(transcript_lineF)
    {
	remove(transcript_line.c_str());
    }
    
    std::string transcript_line_input=Asupp+"transcript_line_input.txt";
    std::ifstream transcript_line_inputF(transcript_line_input.c_str());    
    if(transcript_line_inputF)
    {
	remove(transcript_line_input.c_str());
    }
    
    for (unsigned i = 0; i < boot; i++)
        {       
	    std::stringstream rooted_tree;
	    rooted_tree << Asupp << "rooted_tree_" << i;
	    std::string rooted_treeS = rooted_tree.str();
	    std::ifstream rooted_treeF(rooted_treeS.c_str());    
	    if(rooted_treeF)
	    {
		remove(rooted_treeS.c_str());
	    }
     
	    std::stringstream msa_fasta;
	    msa_fasta << perturbed_msa << "msa_" << i << ".fasta";
	    std::string msa_fastaS = msa_fasta.str();
	    std::ifstream msa_fastaF(msa_fastaS.c_str());    
	    if(msa_fastaF)
	    {
		remove(msa_fastaS.c_str());
	    }
	    
	    std::stringstream msa;
	    msa<< format_msa << "msa_" << i;
	    std::string msaS = msa.str();
	    std::ifstream msaF(msaS.c_str());
	    if(msaF)
	    {
		remove(msaS.c_str());
	    } 
	}

    std::string base_format_file=format_msa+"base_format";
    std::ifstream base_format_fileF(base_format_file.c_str());    
    if(base_format_fileF)
    {
	remove(base_format_file.c_str());
    }
    
    remove(Asupp.c_str());
    remove(perturbed_msa.c_str());
    remove(format_msa.c_str());

    remove(TmpDir.c_str());
}



/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function to create the filtered output file
// Args: name of the output, sequences to analyze, path to the .scores output file, path to the input file, working directory
void filtered_file_distances(const std::string& output_name,std::string& ScoreFilePath, std::string& infilePath,const std::string& WorkDirectory,const std::string& TmpDir,std::map< std::string, int>& CorrNameNum, std::map<int, std::string >& CorrNumName,std::map<std::string, std::vector<int> >& MapGeneTrans,std::vector<int>& TranscriptsList)
{
    ifstream Infile(infilePath.c_str());
    ifstream ScoreFile(ScoreFilePath.c_str());
    
    if (!Infile)
    {
	cerr << " Error in infile: " << infilePath <<endl;
	supp_fichier_distances(TmpDir,output_name);
	exit(37);
    }
    else
    {
	if(!ScoreFile)
	{
	    cerr << " Error in infile: " << ScoreFilePath <<endl;
	    supp_fichier_distances(TmpDir,output_name);
	    exit(38);
	}
	else
	{
	    string FilteredFileName = WorkDirectory + "/"+ output_name + "_filtered.fasta";
	    ofstream FilteredFile(FilteredFileName.c_str(), ios::out | ios::trunc); 
	    if(!FilteredFile)
	    {
		cerr << "cant not create: " << FilteredFileName << endl;
	    }
	    else
	    {	
		// creation of sequence name/score
		string LineNameScore;
		std::map<std::string,float> map_seq_score;
		getline(ScoreFile, LineNameScore);
		while(getline(ScoreFile, LineNameScore))
		    {
			stringstream lineBis(LineNameScore); // Turn the string into a stream.
			std::string SeqName;
			getline(lineBis, SeqName, '\t'); // get sequence name (>ENSP00000375710|Homo_sapiens,>ENSP00000375706|Homo_sapiens,ect)
			std::string score;
			getline(lineBis, score, '\t'); // get sequence score (0.98,0.95,...)
			
			map_seq_score[SeqName]=atof(score.c_str());
		    }
		vector<int> transcriptsKept;
		for (std::map<std::string, vector<int> >::iterator ii = MapGeneTrans.begin(); ii != MapGeneTrans.end(); ii++) 
		    {
			float dist_min=0.0;
			int best_seq;
			std::string best_seq_name;
			for( vector<int>::const_iterator eptr=ii->second.begin(); eptr!=ii->second.end(); eptr++)
			{
			    int SeqNumD=*eptr;
			    std::string SeqNameD;
			    std::map<int, std::string >::iterator itSeqNum;
			    itSeqNum = CorrNumName.find(SeqNumD);
			    if (itSeqNum == CorrNumName.end())
			    {
				cerr << "Error in the CorrNumName key iiii" << SeqNumD << endl;
			    }
			    else
			    {
				SeqNameD=itSeqNum->second;
			    }
			    
			    float dist_seq=map_seq_score.find(SeqNameD)->second;
			    if(dist_min==0.0)
			    {
				best_seq=*eptr;
				dist_min=dist_seq;
			    }
			    else
				if (dist_seq<dist_min)
				{
				    best_seq=*eptr;
				    dist_min=dist_seq;
				}
			}

			transcriptsKept.push_back(best_seq);
		    }
		string FastSeq;
		getline(Infile,FastSeq,'>');
		while(getline(Infile,FastSeq,'>'))
		    {
			stringstream FastSeqSS(FastSeq);
			string SeqNameInfile;
			getline(FastSeqSS,SeqNameInfile);
			SeqNameInfile=">"+SeqNameInfile;
			int SeqNumInfile;
			
			std::map< std::string, int>::iterator itSeqName;
			itSeqName = CorrNameNum.find(SeqNameInfile);
			if (itSeqName == CorrNameNum.end())
			{
			    cerr << "Error in the CorrNameNum key iiii" << SeqNameInfile << endl;
			}
			else
			{
			    SeqNumInfile=itSeqName->second;
			}

			if( (std::find(TranscriptsList.begin(),TranscriptsList.end(),SeqNumInfile)==TranscriptsList.end()) || (std::find(transcriptsKept.begin(),transcriptsKept.end(),SeqNumInfile)!=transcriptsKept.end()) )
			    {
				FilteredFile << ">" << FastSeq << endl;
			    }
		    }
	    }
	    ScoreFile.close();
	    FilteredFile.close();
	}
    }
}



/****************************************************************************************************************************/
/****************************************************************************************************************************/

// Function to remove temporary files
// Args: temporary directory path, bootstrap number and output name
void supp_fichier_distances(const std::string& TmpDir,const std::string& output_name)
{
    std::string Asupp=TmpDir+"/Asupp/";
    std::string perturbed_msa=TmpDir+"/perturbed_msa/";
    std::string format_msa=TmpDir+"/format_msa/";

    std::string InfileS=Asupp+"Infile_number.fasta";
    std::ifstream InfileF(InfileS.c_str());    
    if(InfileF)
    {
	remove(InfileS.c_str());
    }
    
    std::string Name=Asupp+"Name_correspondance.txt";
    std::ifstream NameF(Name.c_str());    
    if(NameF)
    {
	remove(Name.c_str());
    }
   
    std::string align=Asupp+"align_name";
    std::ifstream alignF(align.c_str());    
    if(alignF)
    {
	remove(align.c_str());
    }
    
    std::string alignWithNumber=Asupp+"align_with_nb.fasta";
    std::ifstream alignWithNumberF(alignWithNumber.c_str());    
    if(alignWithNumberF)
    {
	remove(alignWithNumber.c_str());
    }
    
    std::string alignWithNumberGapX=Asupp+"align_with_nb_gap_X.fasta";
    std::ifstream alignWithNumberGapXF(alignWithNumberGapX.c_str());    
    if(alignWithNumberGapXF)
    {
	remove(alignWithNumberGapX.c_str());
    }
    std::string lines_names=Asupp+"lines_names";  
    std::ifstream lines_namesF(lines_names.c_str());    
    if(lines_namesF)
    {
	remove(lines_names.c_str());
    }
      
    std::string matrices=Asupp+"tree.dst";
    std::ifstream matricesF(matrices.c_str());    
    if(matricesF)
    {
	remove(matrices.c_str());
    }
    
    std::string trees=Asupp+"tree.phb";
    std::ifstream treesF(trees.c_str());    
    if(treesF)
    {
	remove(trees.c_str());
    }
    
    std::string matrice_distances=Asupp+"matrice_DistancesScores";
    std::ifstream matrice_distancesF(matrice_distances.c_str());    
    if(matrice_distancesF)
    {
	remove(matrice_distances.c_str());
    }
    
    std::string infile_align=Asupp+"infile_align_desordered.fasta";
    std::ifstream infile_alignF(infile_align.c_str());    
    if(infile_alignF)
    {
	remove(infile_align.c_str());
    }
    
    std::string transcript_line=Asupp+"transcript_line.txt";
    std::ifstream transcript_lineF(transcript_line.c_str());    
    if(transcript_lineF)
    {
	remove(transcript_line.c_str());
    }
    
    std::string transcript_line_input=Asupp+"transcript_line_input.txt";
    std::ifstream transcript_line_inputF(transcript_line_input.c_str());    
    if(transcript_line_inputF)
    {
	remove(transcript_line_input.c_str());
    }

    std::string base_format_file=format_msa+"base_format";
    std::ifstream base_format_fileF(base_format_file.c_str());    
    if(base_format_fileF)
    {
	remove(base_format_file.c_str());
    }
    
    remove(Asupp.c_str());
    remove(perturbed_msa.c_str());
    remove(format_msa.c_str());

    remove(TmpDir.c_str());
}

