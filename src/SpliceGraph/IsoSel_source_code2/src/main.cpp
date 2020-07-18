/*##########################################
#                                          #
#          IsoSel program                  #
#                                          #
#          Héloïse Philippon               #
#          Alexia Souvane                  #
#          Céline Brochier-Armanet         #
#          Guy Perrière                    #
#                                          #
#          LBBE - Lyon - France            #
#          Guy.Perriere@univ-lyon1.fr      #
#                                          #
#          UNIX version, written in C++    #
#                                          #
###########################################*/


#include "fonction.h"
#include "memory.h"

#include <libgen.h>
#include <unistd.h>
#include <stdio.h>
#include <sys/stat.h>

using namespace std;

std::string currentDateTime()
{
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
    return buf;
}
    bool DistancesScores=false;
    bool GapPenalty=false;
    bool ShortPenalty=false;
    
    
int main(int argc, const char* argv[])
{
/*******************************************/
// Variable declaration and initialization //
/*******************************************/

    int nb_thread=1;
    unsigned int nb_boot_int(30);
    string model="LG";
    int miss_dist_est=2;
    std::string output_name = "output" ;
    std::string TranscriptFilePathName = "";
    bool altern = false;
    bool Quiet=false;
    bool only_uniq_prot=true;
    bool autom=false;
    int SeedBoot=0;
    std::string TmpDirUser;
    
    #ifdef _OPENMP
        omp_set_schedule(omp_sched_guided,1);
    #endif
    std::string align_prog="mafft";
    string pathToBatfinder=getenv("ISOSEL_BIN");
    
    // Argument processing
    if (argc < 2)
    {
    cout << "       ********************************     " << endl ;
    cout << "       ****   IsoSel version 1.0   ****     " << endl ;
    cout << "       ********************************     " << endl << endl;
    
        cout << "Usage: run_isosel <infile> [-f <isoforms_locus_file>] [-a alignment_algorithm] [-b <number_of_bootstrap>] [-t <threads_number>] [-m <substitution_model>] [-e <missing_distance_estimation_method>] [-n <output_name>] [-outdir <output_directory>] [-s <number>] [-gap] [-short] [-DS] [-WOT] [-auto] [-quiet]" << endl << endl;
    
    cout << "Basic usage:" << endl;
    cout << "         run_isosel example.fasta " << endl << endl ;
    
    cout << "Input files:" << endl;
    cout << "    <infile>   Multiple sequences input file in fasta format." << endl;
    cout << "    -f         Tab delimited file with sequence names and corresponding coding gene identifier. " << endl << endl;
    
    cout << "Algorithm options:" << endl;
    cout << "    -a         Alignment algorithm. Must be clustalo, maff or muscle (default: mafft). " << endl;
    cout << "    -b         Number of bootstrap replicates. Must be a positive integer (default: 30). " << endl;
    cout << "    -t         Number of threads. Must be a positive integer (default: 1). " << endl;
    cout << "    -m         Substitution model for distance matrices computation. Must be: Blo, Gap, JTT, JTTfast, Kim, LG, Obs, PAM, Pois or WAG (default: LG). " << endl;
    cout << "    -e         Missing distances estimation method. Must be none, MinMax or additive (default: additive) " << endl;
    cout << "    -s         Number used as seed for bootstrap replicates (default: 0) " << endl;
    cout << "    -gap       Option to penalize isoforms introducing gaps in alignment (default: off). " << endl;
    cout << "    -short     Option to penalize isoforms well aligned but too short (default: off). " << endl;
    cout << "    -DS        Option to use the mean of distances as scores (default: off). " << endl;
    cout << "    -WOT       Option to compute scores with all input dataset sequences when the option -DS is used (default: off). " << endl;
    cout << "    -auto      Option to automatically select the best option among -DS, -gap, -short or default parameters (default: off)." << endl;
    cout << "    -quiet     Do not report progress (default: off). " << endl << endl;
    
    
    cout << "Output files:" << endl;
    cout << "    -outdir   Path to the directory results (default: directory of the input file) " << endl;
    cout << "    -n        Output name (default: output) " << endl << endl;
    }
    else
    {
    // Working directory recuperation
        std::string InfilePath = argv[1];
        std::string WorkDir;
        int LastSlash = InfilePath.rfind('/');
        if ( LastSlash == -1 )
        {
            WorkDir="./";
        }
        else
        {
            WorkDir = InfilePath.substr(0, LastSlash)+"/";
        }
        
    // infile opening
            std::ifstream Infile(InfilePath.c_str());

            if (!Infile)
            {
                cout << "Cannot open: " << InfilePath << endl;
                return 0;
            }
         Infile.close();
     
        // other options 
        for (unsigned i = 2; i < argc; i++)
        {       
            if (strstr(argv[i], "-f") != NULL) // If a transcript file is provided
            {
                string AltTranscFilePath = argv[i+1];
                ifstream AltTranscFile(AltTranscFilePath.c_str());
                if (!AltTranscFile)
                {
                    cerr << "Error: unable to open transcript file" << endl;
                    return 2;
                }
                else
                {
                    TranscriptFilePathName = argv[i+1];
                    altern = true;
                }
            }

            else if (strstr(argv[i], "-a") != NULL) // Alignment program choice
            {
		if ( argv[i] == std::string("-a"))
		{
		    if ( argv[i+1] == std::string("muscle") || argv[i+1] == std::string("mafft") || argv[i+1] == std::string("clustalo") )
		    {
			align_prog=argv[i+1];
		    }
		    else
		    {
			cerr << "Error: invalid alignment program. Must be 'muscle' (default) , 'mafft' or 'clustalo'" << endl;
			exit(0);
		    }
		}
		else if(argv[i] == std::string("-auto"))
		{
		    autom = true;
		}
		else{
		    cerr << "Error: invalid option: " << argv[i] << endl;
		    exit(0);
		}
	    }

            else if (strstr(argv[i], "-t") != NULL) // Number of thread processing
            {
        #ifdef _OPENMP
            nb_thread = atoi(argv[i+1]);
            // Checking the the given number
            check_nb_thread(nb_thread);
	    omp_set_num_threads(nb_thread);
            #else 
            fprintf(stderr, "Warning: multi-treading not avalaible, no OpenMP library found\n");
        #endif
            }

            else if (strstr(argv[i], "-b") != NULL) // Number of bootstrap setting
            {
		std::string nb_boot_str= argv[i+1] ;
		if ( (argv[i+1][0]) == '-')
		{
		    cout<< "Error: incorrect number of bootstrap replicates " << endl;
		    exit(1);
		}
		else
		{
		    nb_boot_int = atoi(nb_boot_str.c_str());
		    // Checking the given number
		    check_nb_boot(nb_boot_int);
		}
            }

            
            else if (strstr(argv[i], "-n") != NULL)  // Output Name
            {
		output_name = argv[i+1];
            }
                        
            else if (strstr(argv[i], "-m") != NULL)  // Substitution model
            {
		if ( argv[i+1] == std::string("Blo") || argv[i+1] == std::string("JTT")  || argv[i+1] == std::string("Kim") || argv[i+1] == std::string("LG") || argv[i+1] == std::string("Obs") || argv[i+1] == std::string("Gap") || argv[i+1] == std::string("PAM") || argv[i+1] == std::string("Pois") || argv[i+1] == std::string("WAG") )
		{
		    model= argv[i+1];
		}
		else if(argv[i+1] == std::string("JTTfast") )
		{
		    model= "Gam -a 2.4";
		}
		else
		{
		    cerr << "Error: invalid substitution model. Must be 'Blo', 'Gap', 'JTT', 'JTTfast', 'Kim', 'LG', 'Obs', 'PAM', 'Pois' or 'WAG' (default: LG)." << endl;
		    exit(2);
		}
            }
            
            else if (strstr(argv[i], "-e") != NULL)  // Substitution model
            {
		if (argv[i+1] == std::string("none") )
		{
		    miss_dist_est=0;
		}
		else if (argv[i+1] == std::string("MinMax"))
		{
		    miss_dist_est=1;
		}
		else if (argv[i+1] == std::string("additive") )
		{
		    miss_dist_est=2;
		}
		else
		{
		    cerr << "Error: missing distances estimation method. Must be 'none', 'MinMax' or 'additive' (default: additive)." << endl;
		    exit(3);
		}
            }
            
            
            else if (strstr(argv[i], "-outdir") != NULL) // Output directory
            {
		std::string OutDirPath = argv[i+1];
		int LastSlash = OutDirPath.rfind('/');
		int FirstSlash = OutDirPath.find('/');
		if ( LastSlash == -1 ) // no slash in outdir, example : "Result" or "Result/"
		{
		    WorkDir="./"+OutDirPath+"/";
		}
		else
		    if ( (FirstSlash==1) ) // "./Result" or "./Result/"
		    {
			WorkDir=OutDirPath+"/";
		    }
		    else // "/home/Documents/IsoSel/data/Result" or "/home/Documents/IsoSel/data/Result/"
		    {
			WorkDir = OutDirPath+"/";
		    }

		if (access(argv[i+1], R_OK))
		{
		    cerr << "No output directory: " << WorkDir << endl;
		    exit(4);
		}
	    }
	    else if (strstr(argv[i], "-s") != NULL)
	    {
		if (strstr(argv[i], "-short") != NULL)
		{
		    ShortPenalty = true; 
		}
		else 
		{
		    SeedBoot = atoi(argv[i+1]);
		}
	    }

            else if (strstr(argv[i], "-gap") != NULL)  // Gap Penalty fixed at yes
            {
		GapPenalty = true;
            }
                        
            else if (strstr(argv[i], "-quiet") != NULL)  // Quiet option          
            {
		Quiet = true;
            }
            else if (strstr(argv[i], "-DS") != NULL)  // DS option          
            {
		DistancesScores= true;
            }
            
            else if (strstr(argv[i], "-WOT") != NULL)  // DS plus WOT options          
            {
		if (!altern)
		{
		    cerr << "Transcript file required for the -WOT option" << endl;
		    exit(5);
		}
		
		if (!DistancesScores)
		{
		    cerr << "The -WOT option must be used together with the -DS one" << endl;
		    exit(6);
		}
		only_uniq_prot= false;
            }
        }
      
    #ifdef _OPENMP
	omp_set_num_threads(nb_thread);
    #endif
    
    if (!Quiet)
    {
        cout << "Number of thread(s) used: " << nb_thread << endl; 
    }



//*********************************************************************************************//
//                                Tmp directory setting                                        //
//*********************************************************************************************//
    string TmpDir;
    time_t now= time(NULL);
    std::stringstream now_ss ;
    now_ss << now;
    std::string the_date = now_ss.str();

    if(!access("/dev/shm/", R_OK))
    {
	TmpDir="/dev/shm/"+output_name+"_"+the_date;
    }
    else
    if(!access("/tmp/", R_OK))
    {
        TmpDir="/tmp/"+output_name+"_"+the_date;
    }
    else
    {
        TmpDir=WorkDir+output_name;
    }

//      std::cout << "   TmpDir  =  " << TmpDir << endl;
//     mkdir(WorkDir.c_str(),S_IWRITE);
   
    
//*********************************************************************************************//
//                                    Log file edition                                         //
//*********************************************************************************************//

    string fichier_res_str = WorkDir+output_name+".log";
    ofstream ResFile(fichier_res_str.c_str(), ios::out | ios::trunc);  //openning of the log file
    if(ResFile)  
    {
    
    
//***************************************************************************************************************//
//                                        Beggining of the core program                                          //     
//***************************************************************************************************************//

//*******************************************************************//
//               Tempory folders creating                            //
//******************************************************************//

    mkdir(TmpDir.c_str(),(S_IRUSR | S_IWUSR | S_IXUSR));
    
    string cmd1 = TmpDir+"/Asupp/";
    mkdir(cmd1.c_str(),(S_IRUSR | S_IWUSR | S_IXUSR));
     
    string cmd2 = TmpDir+"/format_msa/";
    mkdir(cmd2.c_str(),(S_IRUSR | S_IWUSR | S_IXUSR));
    
    string cmd3 = TmpDir+"/perturbed_msa/";
    mkdir(cmd3.c_str(),(S_IRUSR | S_IWUSR | S_IXUSR));

    
//*******************************************************************//
//                       Infiles Checking                            //
//*******************************************************************//    
    unsigned nb_seq(0); // number of input sequences
    nb_seq=Seq_name_to_Number(InfilePath,TmpDir,WorkDir,nb_boot_int,output_name,DistancesScores); // sequence names become number

    if (altern)
    {
	checkTranscriptFile(TranscriptFilePathName, InfilePath, TmpDir,nb_boot_int,output_name);
    }
    
   
    if (!Quiet)
    {
	cout << "Number of sequences: " << nb_seq <<endl;
	if (altern)
	{
//             int nb_alt_iso=Isoforms_couting(TranscriptFilePathName,TmpDir,nb_boot_int,output_name);
	    cout << "Number of alternative isoforms: " << Isoforms_couting(TranscriptFilePathName,TmpDir,nb_boot_int,output_name) << endl;
	}
    }
    
    
//*******************************************************************//
//               Reference alignment computing                       //
//*******************************************************************//
    if (!Quiet)
        {
	    cout << "Alignment of " << InfilePath << " in progress" << endl;
        }

    if (align_prog=="muscle")
	{
	    std::stringstream muscleCMD_str;
	    muscleCMD_str << pathToBatfinder << "/" << align_prog << " -in " << TmpDir << "/Asupp/Infile_number.fasta -out " << TmpDir << "/Asupp/infile_align_desordered.fasta -maxiters 2 -quiet";
	    string muscleCMD = muscleCMD_str.str();
	    int output_muscle_reference;
	    output_muscle_reference=system(muscleCMD.c_str());
	    if (output_muscle_reference!=0)
		{
		    cerr << "Error in the alignment of: " << InfilePath << endl;
		    supp_fichier(TmpDir,nb_boot_int,output_name);
		    exit(7);
		}
	    }

    else if (align_prog=="mafft")
        {
	    std::stringstream mafftCMD_str;   
	    mafftCMD_str << pathToBatfinder << "/" << align_prog << " --auto --maxiterate 2 --quiet --thread " << nb_thread << " " << TmpDir << "/Asupp/Infile_number.fasta  > " << TmpDir << "/Asupp/infile_align_desordered.fasta";
	    string mafftCMD = mafftCMD_str.str();
	    int output_mafft_reference;
	    output_mafft_reference=system(mafftCMD.c_str());
	    if (output_mafft_reference!=0)
		{
		    cerr << "Error in the alignment of: " << InfilePath << endl;
		    supp_fichier(TmpDir,nb_boot_int,output_name);
		    exit(8);
		}
        }
    else if (align_prog=="clustalo")
        {
	    std::stringstream clustaloCMD_str; 
	    clustaloCMD_str << pathToBatfinder << "/" << align_prog << " --in=" << TmpDir << "/Asupp/Infile_number.fasta -o " << TmpDir << "/Asupp/infile_align_desordered.fasta --iterations=2 --force --threads=" << nb_thread;
	    string clustaloCMD = clustaloCMD_str.str();
	    int output_clustalo_reference;
	    output_clustalo_reference=system(clustaloCMD.c_str());
	    if (output_clustalo_reference!=0)
		{
		    cerr << "Error in the alignment of: " << InfilePath << endl;
		    supp_fichier(TmpDir,nb_boot_int,output_name);
		    exit(7);
		}
        }

    
    BaseOrder(TmpDir,WorkDir,output_name); // To have the same sequences order in reference MSA as the input dataset and to create the file "/format_msa/+output_name+_base_format"

    
    Seq_num_to_seq_names(TmpDir,WorkDir,output_name,nb_boot_int); // Function for renaming sequences from number to initial name, create the output file "output.aln"
    if (!Quiet)
	{
	    cout << "End of input dataset alignment" << endl;
	}
   

//*******************************************************************//
//                      Variable definitions                         //
//*******************************************************************//

    float** tab_num_base = NULL; //  Array of reference alignment: 0 if a gap, number of residue otherwise
    float* tab_nb_non_gap = NULL; // Array containing the number of residue present at each position of the reference alignment
    
    unsigned nb_site_base(0); // number of sites of the reference alignment



//*******************************************************************//
//                      Input file processing                        //
//*******************************************************************//

    string out_base_file = TmpDir+"/format_msa/base_format";
    ifstream base_file(out_base_file.c_str());
    
    if (!base_file)
        {
            cerr << "No file: " << out_base_file << endl;
	    supp_fichier(TmpDir,nb_boot_int,output_name);
	    exit(9);
        }
    else
        {
	    std::vector<string> align_tab = create_align_tab(base_file);
            nb_seq = align_tab.size();
	    nb_site_base = align_tab[1].size();
	    
	    if (!Quiet)
		{
		    cout << "Number of sites in the reference alignment: " << nb_site_base <<endl;
		}
		
	     // Memory allocation
	    tab_num_base = allocateTwoDim_float(nb_seq,nb_site_base);
	    if (tab_num_base == NULL)
		{
		    cerr << "Memory allocation error" << endl;
		    supp_fichier(TmpDir,nb_boot_int,output_name);
		    exit(10);
		}
            create_tab_num_base(tab_num_base, align_tab, nb_seq, nb_site_base); // Array's filling
	    
	    if(autom)
	    {    
		if (!Quiet)
		{
		    cout << " Automatic selection of IsoSel options... " << endl;
		}
		GapPenalty=false;
		ShortPenalty=false;
		DistancesScores=false;
		if( (nb_seq > 600) || (nb_site_base > 10000) )
		{
		    DistancesScores=true;
                    nb_seq=Seq_name_to_Number(InfilePath,TmpDir,WorkDir,nb_boot_int,output_name,DistancesScores);
		    if (!Quiet)
		    {
			cout << "   Use of the -DS option" << endl;
		    }
		}
		else{
		   if(stats_alignment(tab_num_base, nb_seq, nb_site_base))
		   {
		       GapPenalty=true;
                       if (!Quiet)
                       {
                           cout << "   Use of the -gap option" << endl;
                       }
		   }
		   else{
                       if (!Quiet)
                       {
                           cout << "   Use of the default parameters" << endl;
                       }
                   }
		}
	    }
	    
            tab_nb_non_gap = new float[nb_site_base];
            create_tab_nb_non_gap(tab_nb_non_gap,tab_num_base, nb_seq, nb_site_base); // Array's filling
            base_file.close();
        }

   

    if (DistancesScores==false)
	{
		//**************************************************************//
		//                     Bootstrapping Step                       //
		//**************************************************************//

	    // Formatting for fastdist
	    string base_MSA = WorkDir+output_name+".aln";
	    ifstream base_MSA_file(base_MSA.c_str());
	    rename_for_fastdist(base_MSA_file,TmpDir,DistancesScores);  // file tmpdir/Asupp/align_with_nb.fasta
       
	    // Fastdist for boostrapping and matrice of distance computation
	    ResFile << "Alignment program used: " << align_prog << endl << endl;
    
	    if (!Quiet)
	    {
		cout << "Beginning of distances calculation" << endl;
	    }
       
	    ResFile << "Beginning of distances calculation" << endl;
	    ResFile << "date = " << currentDateTime() << endl << endl;
	    std::stringstream fastdist_cmd_str;
	    fastdist_cmd_str << pathToBatfinder << "/" << "fastdist " << TmpDir << "/Asupp/align_with_nb.fasta -m " << model << " -b " << nb_boot_int << " -t " << nb_thread << " -e " << miss_dist_est << " -s " << SeedBoot << " -n " << TmpDir << "/Asupp/tree -x"; // donne les fichiers ../Asupp/rooted_tree
	    string fastdist_cmd = fastdist_cmd_str.str();
	    int output_fastdist;
	    output_fastdist=system(fastdist_cmd.c_str());
	    if (output_fastdist!=0)
		{
		    cerr << "Error in distance matrices computing" << endl;
		    supp_fichier(TmpDir,nb_boot_int,output_name);
		    exit(11);
		}

	    ResFile << "End of distances calculation and guide trees inference" << endl;
	    ResFile<< "date = " << currentDateTime() << endl << endl ;

	    if (!Quiet)
		{
		    cout << "Beginning of realignement part..." << endl;
		}
		

	    #ifdef _OPENMP
		#pragma omp parallel for 
	    #endif

	    for (unsigned int nb_tree=0 ; nb_tree < nb_boot_int ; nb_tree++)
		{
		    if (!Quiet)
			{
			    cout << "Bootstrap replicate: "  << nb_tree << endl;
			}
		    
		    // Matrice files suppression
		    std::stringstream MatriceSS;
		    MatriceSS << TmpDir << "/Asupp/matrice_" << nb_tree << "_sans_nan";
		    std::string Matrice=MatriceSS.str();
		    remove(Matrice.c_str());

		    // Tree rooting
		    std::stringstream SeaviewCMD;
		    SeaviewCMD<< pathToBatfinder << "/" << "seaview -reroot -root_at_center -outnewick " << TmpDir << "/Asupp/rooted_tree_" << nb_tree << " " << TmpDir << "/Asupp/tree_" << nb_tree <<".phb"; 
		    string Seaview_CMD = SeaviewCMD.str();
		    int output_seaview;
		    output_seaview=system(Seaview_CMD.c_str());
		    if (output_seaview!=0)
			{
			    cerr << "Error in rooting tree " << nb_tree << " using seaview" << endl;
			    supp_fichier(TmpDir,nb_boot_int,output_name);
			    exit(13);
			}

		    // Tree files suppression
		    std::stringstream UnrootedTreeSS;
		    UnrootedTreeSS <<TmpDir << "/Asupp/tree_" << nb_tree << ".phb";
		    std::string UnrootedTree=UnrootedTreeSS.str();
		    remove(UnrootedTree.c_str());

		    // Re alignement step
		    if (align_prog=="muscle")
			{
			    std::stringstream MuscleCMD_str;
			    MuscleCMD_str << pathToBatfinder << "/" <<  align_prog <<" -in " << TmpDir <<"/Asupp/align_with_nb.fasta -out " << TmpDir << "/perturbed_msa/msa_" << nb_tree << ".fasta -usetree_nowarn " << TmpDir << "/Asupp/rooted_tree_" << nb_tree << " -quiet -maxiters 3";
			    string MuscleCMD = MuscleCMD_str.str();
			    int output_muscle_perturbed;
			    output_muscle_perturbed=system(MuscleCMD.c_str());
			    if (output_muscle_perturbed!=0)
				{
				    cerr << "Error in the realignment of bootstrap replicate: " << nb_tree << endl;
				    supp_fichier(TmpDir,nb_boot_int,output_name);
				    exit(14);
				}
			}
		
		    else if (align_prog=="mafft")
			{
			    MafftTreeFormating(TmpDir,nb_tree);
		    
			    std::stringstream mafftCMD_str;  
			    mafftCMD_str << pathToBatfinder << "/" <<  align_prog << " --auto --maxiterate 2 --quiet --thread 1 --treein " << TmpDir << "/Asupp/rooted_tree_mafft_" << nb_tree << " "<< TmpDir << "/Asupp/align_with_nb.fasta  > " << TmpDir << "/perturbed_msa/msa_" << nb_tree << ".fasta";
			    string mafftCMD = mafftCMD_str.str();
			    int output_mafft_perturbed;
			    output_mafft_perturbed=system(mafftCMD.c_str());
			    if (output_mafft_perturbed!=0)
				{
				    cerr << "Error in the realignment of bootstrap replicate: " << nb_tree << endl;
				    supp_fichier(TmpDir,nb_boot_int,output_name);
				    exit(15);
				}
			    std::stringstream rooted_mafft;
			    rooted_mafft << TmpDir << "/Asupp/rooted_tree_mafft_" << nb_tree;
			    std::string rooted_mafftS = rooted_mafft.str();
			    std::ifstream rooted_mafftF(rooted_mafftS.c_str());    
			    if(rooted_mafftF)
				{
				    remove(rooted_mafftS.c_str());
				}
			}
		
		    else if (align_prog=="clustalo")
		    {
			std::stringstream clustaloCMD_str;   
			clustaloCMD_str << pathToBatfinder << "/" <<  align_prog << " --in=" << TmpDir << "/Asupp/align_with_nb.fasta -o " << TmpDir << "/perturbed_msa/msa_" << nb_tree << ".fasta --guidetree-in="<<TmpDir << "/Asupp/rooted_tree_" << nb_tree << " --iterations=2 --force --threads=1";
			string clustaloCMD = clustaloCMD_str.str();
			int output_clustalo_perturbed;
			output_clustalo_perturbed=system(clustaloCMD.c_str());
			if (output_clustalo_perturbed!=0)
			    {
				cerr << "Error in the realignment of bootstrap replicate: " << nb_tree << endl;
				supp_fichier(TmpDir,nb_boot_int,output_name);
				exit(16);
			    }
		    }
		    MsaOrder(TmpDir,nb_tree);
		}
	    
	    #ifdef _OPENMP
		#pragma omp barrier
	    #endif
         
	    if (!Quiet)
		{
		    cout << "Guide tree inference and realignement completed" << endl;
		}
	    ResFile << "Guide tree inference and realignement completed" << endl;
	    ResFile<< "date = " << currentDateTime() << endl << endl;

		//***********************************************************************************//
		//                            Scores computing                                       //
		//***********************************************************************************//
       
                
		//**************************************************************//
		//              Number of scores to compute                     //
		//**************************************************************//
    
	    vector<short> sequence_lines; //vector of sequences for which we have to cumpute score (ie sequences of transcript_file.txt)
	    if (altern == false) // if a transcript file has been provided
		{
		    sequence_lines.resize(nb_seq);
		    for (unsigned i(0) ; i < nb_seq ; ++i)
			{
			    sequence_lines[i] = i;
			}
		}
	    else // if not, we compute all scores
		{
		    sequence_lines=create_species_line(TranscriptFilePathName, InfilePath, TmpDir,nb_boot_int,output_name);
		}

	    // number of scores to compute
	    short nb_species = sequence_lines.size(); 
        
	    // Array of array containing scores of each bootstrap
	    float** tab_residu_score_for1msa[nb_boot_int]; 
    
	    // Array containing scores of each bootstrap
	    float* tab_score_par_seq[nb_boot_int]; 

	    if (!Quiet)
	    {
		cout << "Setting up residues score..." << endl;
	    }
	    ResFile << "Setting up residues score..." << endl;
	    ResFile<< "date = " << currentDateTime() << endl << endl;

	    #ifdef _OPENMP
		#pragma omp parallel for
	    #endif
    
	    for (unsigned int i=0 ; i <  nb_boot_int ; ++i)
		{
		    stringstream ss;
			ss << i;
		    string str = ss.str();
		    string cmd2 = TmpDir+"/format_msa/msa_"+str;

		    ifstream pert_file(cmd2.c_str()); // Opening of the perturbed msa

		    if (!pert_file)
			{
			    cerr << "No such file: " << cmd2 << endl;
			    supp_fichier(TmpDir,nb_boot_int,output_name);
			    exit(17);
			}
		    else
			{
			 /*****************************************************************************************/
			    string* align_tab_pert = new string [nb_seq]; // Memory allocation
			    create_align_tab_pert(pert_file, align_tab_pert);
			    unsigned nb_site_pert = align_tab_pert[1].size(); // number of sites in the perturbed MSA

			/*****************************************************************************************/
			    short** tab_num_pert = allocateTwoDim_short(nb_seq, nb_site_pert); //  Memory allocation for the array containing: 0 if a gap, number of residue otherwise
			    create_tab_num_pert(tab_num_pert, align_tab_pert, nb_seq, nb_site_pert);
			    delete [] align_tab_pert; // Memory deallocation

                        /*****************************************************************************************/
			    tab_residu_score_for1msa[i] = allocateTwoDim_float(nb_species, nb_site_base); // Memory allocation for the array containing residue scores for this bootstrap
			    create_residu_score_for1msa(tab_num_base,tab_num_pert, nb_seq, nb_site_base, tab_residu_score_for1msa[i], sequence_lines, tab_nb_non_gap, GapPenalty); // score computing
			    destroyTwoDim_short(tab_num_pert, nb_seq, nb_site_pert); //  Memory deallocation

                        /*****************************************************************************************/
			    tab_score_par_seq[i] = new float[nb_species];  // Memory allocation for the array containing sequences scores for this bootstrap
			    for (unsigned j(0) ; j < nb_species ; ++j)
				{
				    tab_score_par_seq[i][j] = 0;
				}

			    create_tab_score_seq(tab_score_par_seq[i], tab_residu_score_for1msa[i], nb_species, nb_seq, nb_site_base, sequence_lines, tab_nb_non_gap, ShortPenalty);
			    destroyTwoDim_float(tab_residu_score_for1msa[i], nb_species, nb_site_base); // Memory deallocation

			    pert_file.close();
			}
		} 
        
	    #ifdef _OPENMP
		#pragma omp barrier
	    #endif

	    ResFile << "End of scores calculation" << endl;
	    ResFile<< "date = " << currentDateTime() << endl << endl;
       
	    destroyTwoDim_float(tab_num_base, nb_seq, nb_site_base); // Memory deallocation

    //**************************************************************//
    //                Sequence scores computing                     //
    //**************************************************************//

	    if (!Quiet)
		{
		    cout << "Setting up sequence scores..." << endl;
		}
	    float* score_par_seq_final = new float[nb_species](); // Array of score per sequences

	    short nb(0);
	    for(unsigned i(0) ; i < nb_species ; ++i)
		{
		    for(unsigned j(0) ; j < nb_boot_int ; ++j)
			{
			    score_par_seq_final[i] = score_par_seq_final[i] + tab_score_par_seq[j][i];
			}
		    score_par_seq_final[i] = score_par_seq_final[i] / nb_boot_int;
		}
		
		        
    //***************************************************************************************************************//
    //                                             Output file processing                                            //     
    //***************************************************************************************************************//

          // File with score per sequences
        create_seq_score_file_with_lines(score_par_seq_final, InfilePath, sequence_lines, nb_species, output_name, WorkDir, TmpDir);
        delete [] score_par_seq_final; // Déallocation mémoire
        
         // Filtered file processing
        if (altern == true)
	    {
		string scoreFilePath= WorkDir+output_name+".scores";
		filtered_file(output_name,TranscriptFilePathName,scoreFilePath,InfilePath,WorkDir,TmpDir,nb_boot_int); // create output.filtered
	    }
	}
	
    else // Using the average distances as scores
      {
	    // Formatting for fastdist
	    string base_MSA = WorkDir+output_name+".aln";
	    std::ifstream base_MSA_file(base_MSA.c_str());
	    rename_for_fastdist(base_MSA_file,TmpDir,DistancesScores);  // file tmpdir/Asupp/align_with_nb.fasta
       
	    // Fastdist for boostrapping and matrice of distance computation
	    ResFile << "Alignment program used: " << align_prog << endl << endl;
    
	    if (!Quiet)
	    {
		cout << "Beginning of Distance Scores computation" << endl;
	    }
	    
       nb_boot_int=0; // if the -b option was used by mistake
       
	    ResFile << "Beginning of distances calculation" << endl;
	    ResFile << "date = " << currentDateTime() << endl << endl;
	    std::stringstream fastdist_cmd_str; 
	    fastdist_cmd_str << pathToBatfinder << "/" << "fastdist " << TmpDir << "/Asupp/align_with_nb.fasta -m Gap -b 0 -t " << nb_thread << " -e " << miss_dist_est << " -s " << SeedBoot << " -n " << TmpDir << "/Asupp/tree"; // create the tree_X.phb    
	    string fastdist_cmd = fastdist_cmd_str.str();
// 	    std::cout << "fastdist_cmd_str = " << fastdist_cmd << endl;
	    int output_fastdist;
	    output_fastdist=system(fastdist_cmd.c_str());

 	    if (output_fastdist!=0)
		{
		    cerr << "Error in distance matrices computing" << endl;
		    supp_fichier(TmpDir,nb_boot_int,output_name);
		    exit(18);
		}

	    ResFile << "End of distances calculation" << endl;
	    ResFile<< "date = " << currentDateTime() << endl << endl ;
    
	    // For the NA in the distance matrcies
 	    string MatriceFile= TmpDir+"/Asupp/tree.dst";
 	    ifstream matrices_file(MatriceFile.c_str());

	    //*******************************************************************//
	    //                     Different maps initiation                     //
	    //*******************************************************************//   
	
	
	std::map<std::string,int> MapCorNameNum;
	MapCorNameNum=initMapCorNameNum(TmpDir);
        
        std::map<int, std::string > MapCorNumName;
	MapCorNumName=initMapCorNumName(TmpDir);

	std::map<int, float> TransNumSize;
	initTransNumSize(TmpDir,TransNumSize);
	
	
	std::map<std::string, std::vector<int> > MapGeneTrans;
	MapGeneTrans=initMapGeneTrans(TmpDir,TranscriptFilePathName,MapCorNameNum);
    
	std::map<int, std::vector<int> > MapNumAltTrans;
	MapNumAltTrans=initMapNumAltTrans(MapGeneTrans);
    
	std::vector<int> TransNumList;
	TransNumList=ListTranscriptNumber(MapNumAltTrans);

	
	
	//******************************************************//
	//                Opening file Checking                 //
	//******************************************************//    

	ResFile << "Beginning of score computation..." << endl;
	ResFile<< "date = " << currentDateTime() << endl << endl;
    
	std::string MatriceOkPath =  TmpDir + "/Asupp/tree.dst";
	std::ifstream MatriceOK(MatriceOkPath.c_str());  // to read the file "/Asupp/matrices_sans_nan"
    
	if (!MatriceOK)
	    {
		cerr << " No file: " << MatriceOkPath << endl;
	    }
	else
	    {
		std::string DistFilePath;
		DistFilePath = WorkDir+output_name+".DistancesScores";
		std::ofstream DistFile(DistFilePath.c_str(), ios::out | ios::trunc);     
	    
		if(!DistFile)
		    {
			cerr << "Can not open the file: " << DistFilePath << endl;
			exit(19);
		    }
		else
		    {
			DistFile << "SEQUENCE_NAME" << "\t" << "DISTANCE_SCORES" << "\t" << "SIZE" << endl;
			std::string line;
			getline(MatriceOK, line); // the line with the sequence number
			std::vector<int> AltTrans;

			
			//**************************************************************//
			//                Somme of distance computation                 //
			//**************************************************************//

			    while(getline(MatriceOK, line)) 
				{
				    float sum_altr_size=0;
				    float TransSize=0;
				    std::stringstream lineSpace(line);
				    std::string Valeur1;
				    getline(lineSpace, Valeur1,' '); // line contains the seq number and the matrice line
				    int SeqNum;
				    SeqNum=atoi(Valeur1.c_str());
				    if(std::find(TransNumList.begin(), TransNumList.end(), SeqNum)!=TransNumList.end()) // it is an alternativ transcript
					{
					    std::map<int, std::vector<int> >::iterator itSeqNum;
					    itSeqNum = MapNumAltTrans.find(SeqNum);
					    if (itSeqNum == MapNumAltTrans.end())
						{
						    cerr << "Error in the MapNumAltTrans key " << SeqNum << endl;
						}
					    else
						{
						    AltTrans=itSeqNum->second; 
						    for (int i=0;i!=AltTrans.size();i++)
							{
							    sum_altr_size+=TransNumSize.find(AltTrans[i])->second;
							    if(AltTrans[i]==SeqNum)
							    {
								TransSize=TransNumSize.find(AltTrans[i])->second; // size (aa) of the transcript
							    }
							}
						}
						
					    std::string dist;
					    float DistSum=0;
					    int dist_num=1;
					    
					    if(only_uniq_prot) // only distances to genes without alternative transcripts are sum up
					    {
					    getline(lineSpace, dist,' '); // to suppress the second space after the sequence num
					    while(getline(lineSpace, dist,' '))  // split the matrice file by row so by sequence
						{
						    if (!(std::find(TransNumList.begin(), TransNumList.end(), dist_num)!=TransNumList.end())) // it is not a sequence originating from the same gene
							{
							    stringstream ssDist;
							    float distFloat;
							    ssDist.str(dist);
							    ssDist >> distFloat;
							    DistSum+=distFloat;
							}
						    dist_num+=1;
						}
						DistSum=DistSum/(nb_seq - TransNumList.size() );
						PrintDistFile(MapCorNumName,DistFile,DistSum,SeqNum,TransSize); // to print results in output file .distances
					    }
					    else // all distances to other genes are sum up
					    {
						getline(lineSpace, dist,' '); // to suppress the second space after the sequence num
						while(getline(lineSpace, dist,' '))  // split the matrice file by row so by sequence
						    {
							if (!(std::find(AltTrans.begin(), AltTrans.end(), dist_num)!=AltTrans.end())) // it is not a sequence originating from the same gene
							    {
							    stringstream ssDist;
								float distFloat;
								ssDist.str(dist);
								ssDist >> distFloat;
								DistSum+=distFloat;
							    }
							dist_num+=1;
						    }
						DistSum=DistSum/(nb_seq - AltTrans.size() );
						PrintDistFile(MapCorNumName,DistFile,DistSum,SeqNum,TransSize); // to print results in output file .distances
					    }
					}
				}
		    } // end DistFile opennned
 	    } // end MatriceOk openned

	if (!Quiet)
	    {
		cout << "End of Distance Scores comptation" << endl;
	    }
              
        
//***************************************************************************************************************//
//                                             Output file processing                                            //     
//***************************************************************************************************************//

         // Filtered file processing
        if (altern == true)
        {
        string scoreFilePath= WorkDir+output_name+".DistancesScores";
	filtered_file_distances(output_name,scoreFilePath,InfilePath,WorkDir,TmpDir,MapCorNameNum,MapCorNumName,MapGeneTrans,TransNumList);// create output.filtered
	}
      } // end of DistanceScores

    
//*************************************//
//               The End               //
//************************************//
    ResFile << "End of the IsoSel run" << endl;
    ResFile<< "date = " << currentDateTime() << endl << endl;

    // Deleting temporary file
    if(DistancesScores==true)
	{
	   supp_fichier_distances(TmpDir,output_name);
	}
    else
	{
	    supp_fichier(TmpDir,nb_boot_int,output_name);
	}
       
       if (!Quiet)
       {
       std::cout << endl;
       std::cout << "Results saved in directory: \"" << WorkDir << "\" " << endl;
       std::cout << endl;
       }
      
     ResFile.close();
    }    
 } // end of the IsoSel algorithm
        return 0;
} // end of the main fonction

