#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include <tuple>
#include <limits>


using namespace std;


struct DP_cell{

    int Sscore; // Substitution (S) score
    int Dscore; // Deletion (D) score
    int Iscore; // Insertion (I) score

};


class GlobalLocal{

    
    private:
        
        // Private variables and functions

        vector<vector<DP_cell> > table;
        ofstream output_file_write;
        
        signed int gap_score;
        signed int mismatch_score;
        signed int opening_score;
    
        string seq_1;
        string seq_2;

        int global;
        int match_score;
        int score;


        // Read input file

        void read_input(string file_name){

            ifstream input_data;

            vector<string> sequences;
            
            input_data.open(file_name);

            if (input_data) { 

                string line;
                int line_index = -1;
                string sequence_segment = "";
                
                while (getline(input_data, line)) { 

                    if(line.empty()) { continue; }

                    if(line[0] == '>') {
                        sequences.push_back(sequence_segment);
                        line_index++;
                    } else {
                        sequences[line_index].append(line);
                    }

                }
            } 

            this->seq_1 = sequences[0];
            this->seq_2 = sequences[1];
            input_data.close();

            return ;
        };


        // Substitution score to fill table
        int getS_score(int index_i , int index_j){

            int prev_Sscore = this->table[index_i - 1][index_j -1].Sscore;
            int prev_Dscore = this->table[index_i - 1][index_j -1].Dscore;
            int prev_Iscore = this->table[index_i - 1][index_j -1].Iscore;

            int score = get_MaxValue(prev_Dscore,prev_Iscore,prev_Sscore);


            if (this->global == 0) {
                int zero = 0;
                if(seq_1.at(index_i-1) == seq_2.at(index_j-1)){
                    return max(score + match_score,zero);
                }
                return max(score + mismatch_score,zero);
            }
            if(seq_1.at(index_i-1) == seq_2.at(index_j-1)){
                return score + match_score;
            }

            return score + mismatch_score;

        }

        // Substitution score while retracing
        char get_Sscore_Previous(int index_i , int index_j){

            int prev_Sscore = this->table[index_i - 1][index_j -1].Sscore;
            int prev_Dscore = this->table[index_i - 1][index_j -1].Dscore;
            int prev_Iscore = this->table[index_i - 1][index_j -1].Iscore;

            std::vector<char> types{'d','i','s'};
            
            vector<int> values{prev_Dscore, prev_Iscore, prev_Sscore};

            char type = get_MaxValueType(values,types);

            return type;
        }

        // Deletion score to fill table

        int getD_score(int index_i , int index_j){

            int prev_Sscore = this->table[index_i - 1][index_j].Sscore;
            int prev_Dscore = this->table[index_i - 1][index_j].Dscore;
            int prev_Iscore = this->table[index_i - 1][index_j].Iscore;

            int score = get_MaxValue(prev_Dscore + gap_score,prev_Iscore + opening_score + gap_score,prev_Sscore + opening_score + gap_score);

            return score;
        }

        // Deletion score while retracing
        char get_Dscore_Previous(int index_i , int index_j){

            int prev_Sscore = this->table[index_i - 1][index_j].Sscore;
            int prev_Dscore = this->table[index_i - 1][index_j].Dscore;
            int prev_Iscore = this->table[index_i - 1][index_j].Iscore;

            vector<char> types{'d','i','s'};
            vector<int> values{prev_Dscore + gap_score, prev_Iscore + opening_score + gap_score, prev_Sscore + opening_score + gap_score};


            char type = get_MaxValueType(values, types);

            return type;
        }
        
        // Insertion score to fill table
        int getI_score(int index_i , int index_j){

            int prev_Sscore = this->table[index_i][index_j - 1].Sscore;
            int prev_Dscore = this->table[index_i][index_j - 1].Dscore;
            int prev_Iscore = this->table[index_i][index_j - 1].Iscore;


            int score = get_MaxValue(prev_Iscore + gap_score, prev_Dscore + opening_score + gap_score, prev_Sscore + opening_score + gap_score);

            return score;
        }

        // Insertion score while retracing
        char get_Iscore_Previous(int index_i , int index_j){

            int prev_Sscore = this->table[index_i][index_j - 1].Sscore;
            int prev_Dscore = this->table[index_i][index_j - 1].Dscore;
            int prev_Iscore = this->table[index_i][index_j - 1].Iscore;

            vector<char> types{'d','i','s'};
            vector<int> values{prev_Dscore + opening_score + gap_score,prev_Iscore + gap_score ,prev_Sscore + opening_score + gap_score};
            
            char type = get_MaxValueType(values,types);

            return type;
        }

        // Affinity gap to get max previous cell

        char retrace_case( int index_i , int index_j ,char current_case){
            
            char selected_type;

            if(current_case == 'd'){
                return get_Dscore_Previous(index_i,index_j);
            }
            if(current_case == 'i'){
                return get_Iscore_Previous(index_i,index_j);
            }

            return get_Sscore_Previous(index_i,index_j);
        }

        // Retrace the optimal path

        void retrace(int index_i , int index_j, vector<char> *aligned_str1 , vector<char> *aligned_str2, char cell_max_type){

            if (this->global == 0) {
                if(get_MaxValue(this->table[index_i][index_j].Dscore,this->table[index_i][index_j].Sscore,this->table[index_i][index_j].Iscore) == 0 ){
                    this->output_file_write << "Starting index : "<< index_i << "  " << index_j << endl ;
                    return;

                }
            }

            if (this->global == 1) {
                if(this->table[index_i][index_j].Dscore == 0 && this->table[index_i][index_j].Sscore == 0 && this->table[index_i][index_j].Iscore == 0 ){
                    this->output_file_write << "Starting index : "<< index_i << "  " << index_j << endl ;
                    return;
                }
            }

            char next_case_type;
            
            if(index_i == this->seq_1.length() && index_j == this->seq_2.length()){
                vector<char> types{'d','i','s'};
                vector<int> values{this->table[index_i][index_j].Dscore,this->table[index_i][index_j].Iscore,this->table[index_i][index_j].Sscore};

                cell_max_type = get_MaxValueType(values,types);
            }

            next_case_type = retrace_case(index_i, index_j , cell_max_type);
            
            if(cell_max_type == 'd'){
                
                aligned_str1->push_back(this->seq_1.at(index_i-1));
                aligned_str2->push_back('-');

                return retrace(index_i -1 ,index_j, aligned_str1,aligned_str2,next_case_type);

            }else if(cell_max_type == 'i'){

                aligned_str1->push_back('-');
                aligned_str2->push_back(this->seq_2.at(index_j-1));

                return retrace(index_i ,index_j - 1, aligned_str1,aligned_str2,next_case_type);
            }else{

                aligned_str1->push_back(this->seq_1.at(index_i-1));
                aligned_str2->push_back(this->seq_2.at(index_j-1));

                return retrace(index_i - 1,index_j - 1, aligned_str1,aligned_str2,next_case_type);
            }

            return;
            
        }

                // Maximum value for the provided three values in forward pass.

        int get_MaxValue(int value_1 , int value_2 ,int value_3){

            if (this->global == 1) {
                return max(max(value_1,value_2),value_3);
            }

            int value_4 = 0;
            return max(max(max(value_1,value_2),value_3),value_4);
        }

        // Maximum value in the table

        tuple<int,int> get_Max_All(){

            int max = 0,end_i = 0,end_j = 0,i = 0,j = 0;

            for(vector<DP_cell> v1 : this->table){


                j = 0;
            
                for(DP_cell t : v1){


                    int temp = get_MaxValue(t.Sscore,t.Dscore,t.Iscore);
            
                    if (temp >= max) {


                        max = temp;
                        end_i = i;
                        end_j = j;
                    }
            
                    j++;       
                }
                i++;               
            }
            return make_tuple(end_i,end_j);
        }

        // Type of the maximum value during retracing
        char get_MaxValueType(vector<int> comparables, vector<char> types){

            int max_index = max_element(comparables.begin(), comparables.end()) - comparables.begin();
    
            return types[max_index];
        }


    public:
        
        // Constructor for the class

        GlobalLocal(string sequence_file, int global, int match_score , int mismatch_score , int gap_score, int opening_score ){
        
            read_input(sequence_file);

            if(global == 0 ) {

                // local

                this->output_file_write.open(sequence_file + "-local.txt",ios::app);

            } else if(global == 1 ) {

                // global

                this->output_file_write.open(sequence_file + "-global.txt",ios::app);


            }

            this->output_file_write << "Sequence 1 : ";
            this->output_file_write << this->seq_1 << endl;
            this->output_file_write << "Length : " <<  this->seq_1.length() << endl;
            this->output_file_write << "Sequence 2 : ";
            this->output_file_write << this->seq_2 << endl;
            this->output_file_write << "Length : " <<  this->seq_2.length() << endl;


            
            vector<vector<DP_cell>> table(seq_1.length() +1, vector<DP_cell>(seq_2.length() +1));
            this->table = table;

            this->match_score = match_score;
            this->mismatch_score = mismatch_score;
            this->gap_score = gap_score;
            this->opening_score = opening_score;
            this->global = global;
        }

        
        // Forward pass to fill the table

        void forwardPass(){

            this->output_file_write << endl<< "Output" << endl;
            
            this->output_file_write <<"**********************"<<endl;

            this->table[0][0].Sscore = 0;

            this->table[0][0].Dscore = 0;
            
            this->table[0][0].Iscore = 0;

            int infinity = numeric_limits<int>::infinity();

            for (int index_i = 1 ; index_i <= this->seq_1.length() ;  index_i++ ){

                if (this->global == 1) {
                    this->table[index_i][0].Sscore = 50000 * -1;
                    this->table[index_i][0].Iscore = 50000 * -1;
                    this->table[index_i][0].Dscore =  opening_score + index_i * gap_score;
                }
                if (this->global == 0) {
                    this->table[index_i][0].Sscore = 0;
                    this->table[index_i][0].Iscore = 0;
                    this->table[index_i][0].Dscore = 0;
                }                
            }

            for (int index_j = 1 ; index_j <= this->seq_2.length() ;  index_j++ ){

                if (this->global == 1) {

                    this->table[0][index_j].Sscore = 50000 * -1;

                    
                    this->table[0][index_j].Iscore = opening_score + index_j * gap_score;
                    
                    
                    this->table[0][index_j].Dscore =  50000 * -1;
                }
                if (this->global == 0) {

                    this->table[0][index_j].Sscore = 0;
                    
                    this->table[0][index_j].Iscore = 0;
                    
                    this->table[0][index_j].Dscore = 0;
                
                }      
            }
            

            for (int index1 = 1; index1 <= this->seq_1.length(); index1++){

                for(int index2 =1 ; index2 <= this->seq_2.length(); index2++){
                
                    this->table[index1][index2].Sscore = getS_score(index1,index2);

                    this->table[index1][index2].Dscore = getD_score(index1,index2);
                    
                    this->table[index1][index2].Iscore = getI_score(index1,index2);
                
                }
            
            }
            return;
        }



        // Retrace the optimal path
        void retracePath(){

            int end_i = 0;
            int end_j = 0;



            if (this->global == 0) {
            
            
                tie(end_i,end_j) = get_Max_All();

            }


            if (this->global == 1) {



                end_i = this->seq_1.length();
                end_j = this->seq_2.length();
            }
           
        //    
        //   cout << "End Index : "<< end_i << " " << end_j <<endl;
           
            this->score = get_MaxValue(this->table[end_i][end_j].Dscore,this->table[end_i][end_j].Iscore,this->table[end_i][end_j].Sscore);
            this->output_file_write << "End Index : "<< end_i << " " << end_j <<endl;

            vector<char> aligned_string1;
            vector<char> aligned_string2;

            char type;

            retrace(end_i, end_j, &aligned_string1, &aligned_string2, type);
            
            this->output_file_write <<  endl << "Alignment:  "<<   endl;
            this->output_file_write <<"*****************"<<endl;


            vector<char> pipe_symbol;

            int total_score = 0;
            int extension_count = 0;

            int open_count = 0;
            int mismatch_count = 0;
            int match_count = 0;

            // cout << "Aligned string 1 : " << endl;

            for(int i = aligned_string1.size() - 1; i>= 0 ; i-- ){
                if(aligned_string1[i] ==aligned_string2[i]){
                    
                    pipe_symbol.push_back('|');
                    total_score += this->match_score;
                    match_count +=1;

                }else if (aligned_string1[i] !=aligned_string2[i] && aligned_string1[i] != '-' && aligned_string2[i] != '-'){
                    
                    pipe_symbol.push_back(' ');
                    total_score += this->mismatch_score;

                    mismatch_count +=1;
                }else {
                    pipe_symbol.push_back(' ');
                }
            }

            for(int i = aligned_string1.size() - 1; i>= 0 ; i-- ){

                if(i < aligned_string1.size() && aligned_string1[i] == '-' && aligned_string1[i+1] != '-' ){
                    total_score = total_score + this->opening_score + this->gap_score ;
                    open_count+=1;
                    extension_count+=1;
                }

                if(i < aligned_string1.size() && aligned_string1[i] == '-' && aligned_string1[i+1] == '-' ){
                    total_score += this->gap_score;
                    extension_count+=1;
                }
            }

            // cout << "Total score : " << total_score << endl;

            for(int i = aligned_string2.size() - 1; i>= 0 ; i-- ){

                if(i < aligned_string2.size() && aligned_string2[i] == '-' && aligned_string2[i+1] != '-' ){
                    total_score = total_score + this->opening_score + this->gap_score ;
                    open_count+=1;
                    extension_count+=1;
                }

                if(i < aligned_string2.size() && aligned_string2[i] == '-' && aligned_string2[i+1] == '-' ){
                    total_score += this->gap_score;
                    extension_count+=1;
                }
            }

            this->output_file_write <<endl;
            int size = aligned_string2.size();

            // Printing alignment


            for(int i = 0; i < size ; i+=60 ){
                this->output_file_write << "string 1 : " ;
                
                for(int j = size -i - 1; ((j>= size -i -60) && (j >=0)) ; j-- ){
                    this->output_file_write << aligned_string1[j];
                }

                this->output_file_write <<  endl;
                this->output_file_write << "           " ;

                for(int j = i; ((j<= i+ 60 - 1) && (j < size)) ; j++ ){
                    this->output_file_write << pipe_symbol[j];
                }

                this->output_file_write <<  endl;
                this->output_file_write << "string 2 : " ;

                for(int j = size -i - 1; ((j>= size -i -60) && (j >=0)) ; j-- ){
                    this->output_file_write << aligned_string2[j];
                }
                
                this->output_file_write <<  endl;
            }

            // cout << "Match: " << match_count;


            // Final scores

            this->output_file_write << endl << endl << "Report" <<  endl;
            this->output_file_write <<"************"<<endl;
            
            this->output_file_write << "Table score:" << this->score <<endl;
            this->output_file_write << "Optimal score: " << total_score <<endl;

            this->output_file_write <<endl;
            this->output_file_write << "Number of: "<<endl;
            this->output_file_write << "Match: " << match_count << " ";
            this->output_file_write << " Mismatch: " << mismatch_count<< " ";
            this->output_file_write << " Opening gaps: " << open_count << " ";
            this->output_file_write << " Total gaps: " << extension_count<< " ";


            this->output_file_write <<endl;
            this->output_file_write << "Identities: " << match_count << "/" << size << " (" << round(match_count*100/size) << "%)" << endl;
            this->output_file_write << "Gaps: " << extension_count << "/" << size << " (" << round(extension_count*100/size) << "%)" << endl;

            this->output_file_write.close();
        }



};


