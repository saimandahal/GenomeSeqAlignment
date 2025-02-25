#include <iostream>
#include <stdlib.h>
#include <fstream>
#include "GlobalLocal.cpp"

using namespace std;


// Change the parameter.config file to change the parameters. During compilation use --- ./exec seq_file_name alignment_choice(0: local, 1: global)>
// Output are saved in seq_file_name.global.out or seq_file_name.local.out


int main(int arguments, char** args_taken)
{


   if (arguments < 3) {

      cout << " Format: <executable_file seq_file_name alignment_choice(0: local, 1: global)>  " << "\n";

   }

   string seq_file_input = args_taken[1];

   int alignment_taken = atoi (args_taken[2]);
   
   string param_file = "parameters.config";

   ofstream output_file_write;
   string align_name;

   if(alignment_taken == 0 ) {


      output_file_write.open(seq_file_input + "-local.txt");
      align_name = "Local";
   
   
   } else if(alignment_taken == 1 ) {
      

      output_file_write.open(seq_file_input + "-global.txt");
      align_name = "Global";

   } else{
      cout << "Invalid input (0: local, 1: global) ";
   }

   // Starting to write in output file

   output_file_write << "Input" <<endl;
   output_file_write <<"**********************"<<endl;
   output_file_write<<endl;


   output_file_write << "Alignment:\t" << align_name << endl;

   int match_score, mismatch_score, gap_score, opening_score;
   
   ifstream input_data;
   input_data.open(param_file);

   string params;
   int score;

   while (input_data >> params >> score ) { 

      if(params == "match"){
         match_score = score;
      }
      if(params == "mismatch"){
         mismatch_score = score;
      }
      if(params == "h"){
         opening_score = score;
      }
      if(params == "g"){
         gap_score = score;
      }

   }

   output_file_write<< "Scores: "<<endl;
   output_file_write << "Match: " << match_score;
   output_file_write << "Mismatch: " << mismatch_score;
   output_file_write << " Gap opening: " << opening_score;
   output_file_write << " Gap: " << gap_score;
   output_file_write << endl;

   output_file_write <<"**********************"<<endl;

  
   input_data.close();

   // cout<<"I am here"<<endl;


   // Alignment code

   if (alignment_taken == 0) {

      // Local alignment
      
      GlobalLocal alignment_g(seq_file_input, alignment_taken, match_score , mismatch_score, gap_score, opening_score);
      alignment_g.forwardPass();
      alignment_g.retracePath();  

   } else if (alignment_taken == 1) {

      // Global alignment

      GlobalLocal alignment_g(seq_file_input, alignment_taken, match_score , mismatch_score, gap_score, opening_score);
      alignment_g.forwardPass();
      alignment_g.retracePath();   

   } else {

      cout << "Invalid (0: local, 1: global) ";

   }

   output_file_write.close();
   
   return 0;

}
