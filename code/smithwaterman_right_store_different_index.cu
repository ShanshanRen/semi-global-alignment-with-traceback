#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <cuda.h>
#include <stdint.h>
#include <math.h>
#include <unistd.h>
#include <omp.h>
#include<limits>
 double diff(timespec start, timespec end)
 {
 	double a=0;
        if((end.tv_nsec-start.tv_nsec)<0)
        {
        a=end.tv_sec-start.tv_sec-1;
        a+=(1000000000+end.tv_nsec-start.tv_nsec)/1000000000.0;
        }
        else
        {
        a=end.tv_sec-start.tv_sec+(end.tv_nsec-start.tv_nsec)/1000000000.0;

        }
        return a;

}

struct NUM_ADD
{
                short2 read_reference_number;
                int address_array;
};


  __global__ void  calculate_cigar( int size, char * data, NUM_ADD *num_add,int * result, char * cigar,int * cigar_int,int * direction) //, char * result
{
	int offset=blockIdx.x;
	__shared__ short2 read_reference_number;
	__shared__ char * read_base_array;
	__shared__ char4 * reference_base_array;
	__shared__ int mismatch;
	__shared__ int  match;
	__shared__ int open;
	__shared__ int extend;
	__shared__ int * direction_index;
	__shared__ char * cigar_store;
	__shared__ int *cigar_int_store;
	while(offset<size)
         {
		if( threadIdx.x==0)
                {
                read_reference_number=num_add[offset].read_reference_number;
                read_base_array=(char *) (data+num_add[offset].address_array);
	        reference_base_array=(char4 *) (read_base_array+(read_reference_number.x+4-1)/4*4);
		direction_index=(int *) (direction+offset*800*1100);  
		cigar_store=(char *) (cigar+offset*sizeof(char)*128);
		cigar_int_store=(int *) (cigar_int+offset*128);
	     }
                __syncthreads();

 		__shared__ char reference_base_in_char[600];
                int hh=(read_reference_number.y+4-1)/4;
                int tt=(hh+blockDim.x-1)/blockDim.x;
                for(int ii=0;ii<tt;ii++)
                {
                        int aa=threadIdx.x+ii*blockDim.x;
                        if(aa< hh)
                        {
                        char4 reference_base_in_thread;
                        reference_base_in_thread=reference_base_array[aa]; //Is it right to get data from global memory
                        reference_base_in_char[aa*4]=reference_base_in_thread.x;
                        reference_base_in_char[aa*4+1]=reference_base_in_thread.y;
                        reference_base_in_char[aa*4+2]=reference_base_in_thread.z;
                        reference_base_in_char[aa*4+3]=reference_base_in_thread.w;
			 }
                }
		
		__shared__ int MM[130]; 
                __shared__ int gap_h[130]; //insertion
                __shared__ int gap_size_h[130];  //insertion
	 	__shared__ int result_col;
		__shared__ int result_row;
		__shared__ int result_col_index;
		__shared__ int result_row_index;
		__shared__ int final_result;
		__shared__ int final_i;
		__shared__ int final_j;
		__shared__ char cigar_m[128];
		__shared__ int cigar_int_m[128];

                if(threadIdx.x==0)
                {
                        MM[0]=0;
                        gap_h[0]=-1000000000;//std::numeric_limits<int>::min()/2;
                        gap_size_h[0]=0;
			match=200;
			mismatch=-150;
			open=-260;
			extend=-11;
			result_col=-1000000000;//std::numeric_limits<int>::min()/2;
			result_row=-1000000000;//std::numeric_limits<int>::min()/2;
                //	for(int i=0;i<read_reference_number.y;i++)
		//	printf("%c",reference_base_in_char[i]);
		//	printf("\n");
		//	for(int i=0;i<read_reference_number.x;i++)
		//	printf("%c",read_base_array[i]);
		}

		 __syncthreads();


		int read_number=read_reference_number.x;
                int round_size;
		int round=(read_reference_number.x+blockDim.x-1)/blockDim.x;
		for(int i=0;i<round;i++)
                {
                        round_size=(read_number>blockDim.x)?blockDim.x: read_number;
                        read_number=(read_number>blockDim.x)?read_number-blockDim.x:0; // read_num is the remaining length at this round
			
			char read_base;
			 if(threadIdx.x<round_size ) // tid is from 0 ~ round_size-1
                        {
                                read_base=read_base_array[threadIdx.x+blockDim.x*i];
			}

			int gap_v=-1000000000;//std::numeric_limits<int>::min()/2;;
			int gap_size_v=0; //Deletion
 			int M=0; //now 
                        int step_right; //now
			int ki=0;//insertion  h  negetive
			 //deletion  v
			int MMM=0; 
			int curBtrack;
                                          
			int current_reference_id=0;
                        for(int j=0;j<round_size+read_reference_number.y-1;j++)
                        {
                                int aa=j-threadIdx.x;
                                if( aa>=0 && (current_reference_id<read_reference_number.y))
                                {
					int prev_gap=M+open; //M which is cacluated by last step in the same thread
					gap_v+=extend;
					if(prev_gap>gap_v)
					{
						gap_v=prev_gap;
						gap_size_v=1;
					}	
					else
						gap_size_v++;
					
					char reference_base_each=reference_base_in_char[current_reference_id];
				
				 	current_reference_id++;
					M=MMM+(read_base==reference_base_each? match:mismatch);
					prev_gap=MM[threadIdx.x]+open;
					step_right=gap_h[threadIdx.x]+extend;
					if(prev_gap>step_right)
					{
						step_right=prev_gap;
						ki=1;
					}	
					else
						ki=gap_size_h[threadIdx.x]+1;

					bool diag=(M>=gap_v)&&(M>=step_right);
					if(diag)
					{
					  curBtrack=0;
					}
					else
					if(step_right>=gap_v)
					{
						 M=step_right;
						 curBtrack=0-ki;
					}
					else
						{
							M=gap_v;
							curBtrack=gap_size_v;
						}
					MMM=MM[threadIdx.x];
					direction_index[800*j+threadIdx.x]=curBtrack;
				
				//if(threadIdx.x==read_reference_number.x-3)
				//printf("%p %d ", &direction_index[800*j+threadIdx.x],curBtrack);
 				
				if(current_reference_id==read_reference_number.y)
				{	
					if(M>=result_row)
					{
						result_row=M;
						result_row_index=threadIdx.x;  //
					}
					//printf("%d %d  %d  %d %d \n",read_reference_number.y,M,result_row,result_row_index,threadIdx.x);
				}
                         	if(threadIdx.x==round_size-1 && i==round-1)
                                {
						if(M>=result_col)
						{
							result_col=M;
							result_col_index=current_reference_id;	// +1					
						}						
				}



			       }
                        	
				__syncthreads(); //to make sure that the former value of MM[threadIdx.x+1] are used by other threads.
                                MM[threadIdx.x+1]=M;
                                gap_h[threadIdx.x+1]=step_right;
                                gap_size_h[threadIdx.x+1]=ki;

                                                               __syncthreads(); // there should be two synthreads(); // to make sure that all of MM[threadIdx.x+1] have get a new value before M,D and I changed.
                        }
                }
		
		__shared__ int state;//0  match;  1 mistmatch; 2 inseriton;  3  deletion
		__shared__ int cigar_index;
		__shared__ int segment_length;
		__shared__ int btr;
               	__shared__ int new_state;
		__shared__ int step_length;
		 if(threadIdx.x==round_size-1)
                {
			//printf("%d %d %d %d\n", result_row,result_col, result_row_index,result_col_index);
                        if(result_row>result_col||result_row==result_col&&(read_reference_number.x-result_row_index-1)>(read_reference_number.y-result_col_index))
			{
				final_result=result_row;			
				final_i=read_reference_number.y-1;
				final_j=result_row_index;
				segment_length=read_reference_number.x-1-result_row_index;
			}
			else
			{
				final_result=result_col;
				final_i=result_col_index-1;
				final_j=read_reference_number.x-1;
				segment_length=0;
			}
			result[offset*3]=final_result;
			//printf("%d\n",final_result);
               		cigar_index=0;	
			if(segment_length>0)
			{
			cigar_m[cigar_index]='S';
			cigar_int_m[cigar_index]=segment_length;
			segment_length=0;
			cigar_index++;
			}
			
			//printf("\n %d %d\n", final_i,final_j);
			state=0;
			do
			{
				btr=direction_index[(final_i+final_j)*800+final_j];
				if(btr>0)
				{
					new_state=3;
					step_length=btr;
				}
				else
				if(btr<0)
				{
					new_state=2;
					step_length=0-btr;		
				}
				else	
				{
					new_state=0;
					step_length=1;
				}
				
		//		printf(" %d %d %d %d\n", state,new_state,final_i,final_j);
				if(new_state==0)
				{	final_i--;
					final_j--;
				}
				else
				if(new_state==2)
					final_j-=step_length;
				else
					final_i-=step_length;
			
				if(state==new_state) segment_length+=step_length;
				else
				{
				if(state==0)  cigar_m[cigar_index]='M';
				if(state==2)  cigar_m[cigar_index]='I';
				if(state==3)  cigar_m[cigar_index]='D';
                       		 cigar_int_m[cigar_index]=segment_length;
                       		 segment_length=step_length;
                       		 cigar_index++;
			 	 state=new_state;
				}
	
			}while(final_i>=0&&final_j>=0);
			
			if(state==0)  cigar_m[cigar_index]='M';
			if(state==2)  cigar_m[cigar_index]='I';
			if(state==3)  cigar_m[cigar_index]='D';
                       	cigar_int_m[cigar_index]=segment_length;
                       	cigar_index++;
			if(final_j>=0) 
			{	
				cigar_m[cigar_index]='S';
				cigar_int_m[cigar_index]=final_j+1;
				cigar_index++;
			}	

			result[offset*3+1]=final_i+1;
			result[offset*3+2]=cigar_index;
	/*		for(int i=cigar_index-1;i>=0;i--)
			{
			printf("%d%c",cigar_int_m[i],cigar_m[i]);
			}
*/
		 }
		 __syncthreads();
		if(threadIdx.x<cigar_index && cigar_index<=blockDim.x)
		{
	//	if(threadIdx.x==0)
	//		printf("%c %d\n",cigar_m[cigar_index-1-threadIdx.x], cigar_int_m[cigar_index-1-threadIdx.x]);
		cigar_store[threadIdx.x]=cigar_m[cigar_index-1-threadIdx.x];
		cigar_int_store[threadIdx.x]=cigar_int_m[cigar_index-1-threadIdx.x];
	//	if(threadIdx.x==0)
	//		printf("%c %d\n", cigar_store[threadIdx.x],cigar_int_store[threadIdx.x]);
		
		}

		offset+=gridDim.x;
	}
}

struct InputData
{
char read_base[600];
char reference_base[600];
};


int main()
{
	FILE * file;
        file=fopen("pairHMM_input.txt","r");
        int size;
        fscanf(file,"%d",&size);
	double  computation_time=0;//total_time=0;
	timespec start,finish;	 
       //printf("size=%d\n",size);
         while(!feof(file))
        {
                InputData * inputdata=(InputData* )malloc(size*(sizeof(InputData)));
                for(int i=0;i<size;i++)
                {
                        fscanf(file,"%s ",inputdata[i].reference_base);
                        fscanf(file,"%s ",inputdata[i].read_base);
                }
                
		//data preparation.
		char * data_h_total=(char*)malloc(size * 600* sizeof (char)*2);
		NUM_ADD * data_num_add=(NUM_ADD *) (data_h_total);
		char * data_h=data_h_total+size*sizeof(NUM_ADD);  // it is 64*x .thus we donot need to worry about align
		int data_size=0;
		char * data_d_total;		
		cudaMalloc( (char **) &data_d_total, size * 600* sizeof (char)*2+sizeof(int)*size*3);
		int * result_h=(int *) malloc(sizeof(int)*size*3);
	 	char * cigar_h=(char *) malloc(sizeof(char)*size*128);
		int * cigar_int_h=(int *) malloc(sizeof(int)*size*128);
		for(int i=0;i<size;i++)
		{

			char4 reference_tep[150];
			int read_len=strlen(inputdata[i].read_base);
			int ref_len=strlen(inputdata[i].reference_base);
			int new_len=(ref_len+4-1)/4;
			for(int j=0;j<new_len;j++)
		        {
		        	reference_tep[j].x=inputdata[i].reference_base[j*4];
		                if(j*4+1<ref_len)
		                reference_tep[j].y=inputdata[i].reference_base[j*4+1];
		                if(j*4+2<ref_len)
		                reference_tep[j].z=inputdata[i].reference_base[j*4+2];
		                if(j*4+3<ref_len)
		                reference_tep[j].w=inputdata[i].reference_base[j*4+3];                   
		         }
		
			data_num_add[i].read_reference_number.x=read_len;
			data_num_add[i].read_reference_number.y=ref_len;
			data_num_add[i].address_array=data_size;

			memcpy(data_h,inputdata[i].read_base,read_len);
			data_h+=(read_len+4-1)/4*4;
			data_size+=(read_len+4-1)/4*4;

			memcpy(data_h,reference_tep,sizeof(char4)* new_len);
		        data_h+=new_len*sizeof(char4);
		        data_size+=new_len*sizeof(char4);
		}
		
		int data_size_to_copy=data_size+size*sizeof(NUM_ADD);

		cudaMemcpy(data_d_total,data_h_total,data_size_to_copy,cudaMemcpyHostToDevice);
		NUM_ADD * num_add_d=(NUM_ADD *) (data_d_total);
		char * data_d=data_d_total+sizeof(NUM_ADD)*size;
		int * result_d=(int *) (data_d_total+data_size_to_copy);
		
		char * cigar;
		cudaMalloc( (char **) &cigar, size * (128* sizeof (char)+128*sizeof(int)));
		
		int * cigar_int=(int *) (cigar+size*128*sizeof(char));
		int * direction;

		cudaMalloc( (int **) & direction, size * (800*1100* sizeof (int)));
	
		dim3 block(128);
		dim3 grid(size);
		clock_gettime(CLOCK_MONOTONIC_RAW,&start);
		calculate_cigar<<<grid,block>>> (size,data_d,num_add_d,result_d,cigar,cigar_int,direction); //result
		cudaMemcpy(result_h,result_d,size*sizeof(int)*3,cudaMemcpyDeviceToHost);
		cudaMemcpy(cigar_h,cigar,128*sizeof(char)*size, cudaMemcpyDeviceToHost);
		cudaMemcpy(cigar_int_h,cigar_int,128*sizeof(int)*size,cudaMemcpyDeviceToHost);		

		clock_gettime(CLOCK_MONOTONIC_RAW,&finish);
                computation_time+=diff(start,finish);
	
		for(int i=0;i<size;i++)
		{
			printf("%d %d\n",result_h[i*3],result_h[i*3+1]);
			printf("[");
			for(int j=0;j<result_h[i*3+2];j++)
			{
			if(j!=0) printf(", ");
			printf("%d%c",cigar_int_h[128*i+j],cigar_h[128*i+j]);
			}
			printf("]\n");
		}
		cudaFree(direction);
		free(data_h_total);
		cudaFree(data_d_total);
		free(inputdata);
		cudaFree(cigar);
		free(cigar_int_h);
		free(cigar_h);
                fscanf(file,"%d",&size);
        }

 	printf(" computation_time= %e  total_time=%e \n",computation_time,0);



        return 0;
}



