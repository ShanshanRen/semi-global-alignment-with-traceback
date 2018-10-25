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

#define MEM_OFFSET gridDim.x*blockDim.x  //the number of threads in the grid
#define BACK(x,y)   back[startPosA[blockThread] + ( ( ((y) + 8) / 8) * rowWidth + (x) + 1 ) * MEM_OFFSET]
#define Y_STEPS 8
#define BLOCK_SIZE   128   //128   //blockDim.x
#define INT_INT -2147483647

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


//calculate_cigar( int size, char * data, NUM_ADD *num_add,int4 * result, int * direction) 
//rowWidth=640

__global__ void Calculate_ben(int size, char *data, NUM_ADD * num_add, short4 *result, int2* AF_maxXY, unsigned int* back,  short rowWidth)
{
	int offset=blockIdx.x*blockDim.x+threadIdx.x;
	short2 lengthXY;
	char * read_base_array;
	char4 * reference_base_array;
	int  mismatch;
	int  match;
	int gapOp;
	int  gapEx;
	__shared__  int startPosA[128];  //blockDim.x
	while(offset<size)
	{
		int blockThread=threadIdx.x;
		//if(offset>=2943) printf("id=%d\n",offset);
		match=200;
		mismatch=-150;
		gapOp=-260;
		gapEx=-11;
		lengthXY=num_add[offset].read_reference_number;
		read_base_array=(char *) (data+num_add[offset].address_array);
		reference_base_array=(char4 *) (read_base_array+(lengthXY.x+127)/128*128);
		// direction_index=(short2 *)(direction+offset*640*1100);  
		startPosA[threadIdx.x] = offset;  // startPosA 是放的是什么？
		//startPosA 在原程序中是代表 thread在grid中的编号
		//	printf("%d %d\n", lengthXY.x, lengthXY.y);

		//下面是初始化

		//initialization of the -1 row in A matrix
		// - 2 bytes for element of A matrix
		// - 2 bytes for element of F matrix
		//还是不知道AF_maxXY 是放什么的？
		for(short x = 0; x < lengthXY.x; x++)
		{
			int2 tmp;
			//(x + 1) because the first element should be -gapEx
			tmp.x = 0;
			tmp.y = INT_INT - gapEx;
			AF_maxXY[startPosA[threadIdx.x] + x * MEM_OFFSET ] = tmp;  //because of this operation, the total number of threads in the grid should be greater than size. Otherwise, AF would re-written.
			//MEM_OFFSET在原程序中是 grid的总的thread的个数。
			//AF_maxXY[startPosA[blockThread] + x * MEM_OFFSET] = tmp;

			//fill the -1 row of "back" array
			BACK(x,-1) = 9; //0000 0000 0000 0000 0000 0000 0000 1001 == 9
		}

	  //	if(offset>=2943) printf("id=%d\n",offset);
		//	printf("%d %d\n", lengthXY.x, lengthXY.y);
		//fill the -1 column of "back" array
		for(short y = 0; y < lengthXY.y; y+=Y_STEPS)
		{
			
	//	if(offset>=2943) printf("id=%d %d  %d\n",offset,y,startPosA[threadIdx.x] + ( ( ((y) + 8) / 8) * rowWidth + (-1) + 1 ) * MEM_OFFSET);
			BACK(-1,y) = 1717986918; //0110 0110 0110 0110 0110 0110 0110 0110 = 1717986918
		}
		BACK(-1,-1) = 0; //stop element

		
		//one element of AE_shared consist of:
		// - one A element
		// - one E element
		__shared__ int2 AE_shared[Y_STEPS][BLOCK_SIZE];
		//elements of Y sequence go to sharedYSeq
		__shared__ char4 sharedYSeq[Y_STEPS/4][BLOCK_SIZE];


		int2 AF_current;
		AF_current.x = 0;

		__shared__ int2 ymin_score[BLOCK_SIZE]; //stores ymin and score
		ymin_score[threadIdx.x].y = 0;

		__shared__ short4 maxXY[BLOCK_SIZE];
		maxXY[threadIdx.x].x = lengthXY.x - 1;
		maxXY[threadIdx.x].y = 0;
		maxXY[threadIdx.x].z = 0;

		// |
		// |
		// |
		// V
		for (short y = 0; y < lengthXY.y; y += Y_STEPS)
		{
			//printf("%d\n",y);
			int2 A_init_upleft;
			A_init_upleft.x = 0;

			//initialialization of the -1 column in A matrix
			// - one element of A matrix
			// - one element of E matrix
			for (short i = 0; i < Y_STEPS; i++)
			{
				int2 tmp;
				tmp.x = 0;
				tmp.y = INT_INT - gapEx;
				AE_shared[i][threadIdx.x] = tmp;
			}


			//we read elements of the Y sequence
			for (short i = 0; i < Y_STEPS/4; i++)
			{
				sharedYSeq[i][threadIdx.x] = reference_base_array[y/4+i];
				//PACK_BYTES(tex1Dfetch(texSeqsY, startY + y + i*4 + 0),
				//                                       tex1Dfetch(texSeqsY, startY + y + i*4 + 1),
				//                                      tex1Dfetch(texSeqsY, startY + y + i*4 + 2),
				//                                     tex1Dfetch(texSeqsY, startY + y + i*4 + 3));
				//printf("%c %c %c %c\n", sharedYSeq[i][threadIdx.x].x,sharedYSeq[i][threadIdx.x].y,sharedYSeq[i][threadIdx.x].z,sharedYSeq[i][threadIdx.x].w);
			}

			ymin_score[threadIdx.x].x = min(Y_STEPS, lengthXY.y - y); //(i < Y_STEPS) && (i + y < lengthY)

			//------>
			for (short x = 0; x < lengthXY.x; x++)
			{
				//actual up_left gets a value of recent read value from the global memory
				//and actual read value is stored in first two bites of A_upleft
				A_init_upleft.y = A_init_upleft.x;

				char2 XYSeq;
				XYSeq.x = read_base_array[x];
				//	if(y==0) printf("%c\n",XYSeq.x);
				//read from global memory
				int2 AF_up = AF_maxXY[startPosA[threadIdx.x] + x * MEM_OFFSET];

				//A_init -> up element read in previous iteration from global memory (up-left)
				A_init_upleft.x = AF_up.x;

				int2 AE_left;
				int E_current;
				int similarity;
				unsigned int back8 = 0;
				//  |  /|  /|
				//  | / | / |
				//  |/  |/  V
				//  |  /|  /|
				//  | / | / |
				//  |/  |/  V
				for(short i = 0; i < ymin_score[threadIdx.x].x; i++)
				{
					AE_left = AE_shared[i][threadIdx.x];


					// XYSeq.y = sharedYSeq[i/4][threadIdx.x].x,y,z,w;
					if(i%4==0)
						XYSeq.y = sharedYSeq[i/4][threadIdx.x].x;
					if(i%4==1)
						XYSeq.y = sharedYSeq[i/4][threadIdx.x].y;
					if(i%4==2)
						XYSeq.y = sharedYSeq[i/4][threadIdx.x].z;
					if(i%4==3)
						XYSeq.y = sharedYSeq[i/4][threadIdx.x].w;
					//(sharedYSeq[i/4][threadIdx.x] >> (((15-i)%4) * 8)) & 0xFF;


					//similarity = substitutionMatrix[XYSeq.y*lettersCount + XYSeq.x];
					similarity =   (XYSeq.x==XYSeq.y? match:mismatch);
					similarity += A_init_upleft.y;

					E_current = max(AE_left.y + gapEx, AE_left.x + gapOp);
					AF_current.y = max(AF_up.y + gapEx, AF_up.x + gapOp);

					AF_current.x = max(E_current, AF_current.y);
					AF_current.x = max(AF_current.x, similarity);

					//"back" array
					back8 <<= 1;
					//back8 |= ((AF_current.x==E_current) && (AF_current.x!=AF_current.y)) || (AF_current.x==similarity); //if go left
					back8 |= (AF_current.x==E_current)  || (AF_current.x==similarity); //if go left
					back8 <<= 1;
					//back8 |= (AF_current.x==AF_current.y) || (AF_current.x==similarity); //if go up
					back8 |=( (AF_current.x==AF_current.y)&& (AF_current.x!=E_current)) || (AF_current.x==similarity); //if go up
					back8 <<= 1;
					back8 |= (AF_current.y == (AF_up.y + gapEx)); //if continue up
					back8 <<= 1;
					back8 |= (E_current == (AE_left.y + gapEx)); //if continue left

					//initialize variables for next iterations
					int2 AE_tmp;
					AE_tmp.x = AF_current.x;
					AE_tmp.y = E_current;
					AE_shared[i][threadIdx.x] = AE_tmp;
					A_init_upleft.y = AE_left.x;
					AF_up = AF_current;
					//	printf("%d ",AF_current.x);
				}  //end of i
				//printf("\n");
				//we want the last row of back8 to be completed
				back8 <<= 4 * (Y_STEPS - ymin_score[threadIdx.x].x);

				//write variables to global memory for next loop
				AF_maxXY[startPosA[threadIdx.x] + x * MEM_OFFSET] = AF_current;
				BACK(x,y) = back8;

				//looking for max element in the last row
				if( (y + ymin_score[threadIdx.x].x) == lengthXY.y )
				{
					if (AF_current.x > ymin_score[threadIdx.x].y)
					{
						maxXY[threadIdx.x].x = x;
						maxXY[threadIdx.x].y = y + ymin_score[threadIdx.x].x - 1; //why minus 1???? Because 0+8=8,it should be 7.
						maxXY[threadIdx.x].z=lengthXY.x-1-x;
					}
					//if result== last row
					//result4.x=read_reference_number.y-1;
					//result4.y=result_row_index;///result_row_index is the threadIdx.x,which is x.
					//result4.z=read_reference_number.x-1-result_row_index;
					ymin_score[threadIdx.x].y = max(ymin_score[threadIdx.x].y, AF_current.x);
				}

			} //end of x

			//looking for max element in the last column
			for(short i = 0; i < ymin_score[threadIdx.x].x; i++)
			{
				if (AE_shared[i][threadIdx.x].x > ymin_score[threadIdx.x].y||AE_shared[i][threadIdx.x].x==ymin_score[threadIdx.x].y&& maxXY[threadIdx.x].z>(lengthXY.y-(y+i)-1))
				{
					maxXY[threadIdx.x].x = lengthXY.x - 1; //
					maxXY[threadIdx.x].y = y + i;
					maxXY[threadIdx.x].z=0;
				}
				//result4.x=result_col_index;   //result_col_index is the y.
				//result4.y=read_reference_number.x-1;
				//result4.z=0;
				ymin_score[threadIdx.x].y = max(ymin_score[threadIdx.x].y, AE_shared[i][threadIdx.x].x);
			}
		}//end of y
		//        maxXY[threadIdx.x].w=ymin_score[threadIdx.x].y;
		//here write result (AF_current) to global memory
		//   scores[startPosA[blockThread]] = ymin_score[blockThread].y;
		// AF_maxXY[startPosA[threadIdx.x]] = maxXY[threadIdx.x];
		result[offset]=maxXY[threadIdx.x];
		//	printf("%d %d %d %d %d\n",offset,result[offset].x,result[offset].y,result[offset].z,result[offset].w);

		offset+=gridDim.x*blockDim.x;
	}


}

#undef BACK
#define BACK(x,y)   back[startPosA + ( ( ((y) + 8) / 8) * rowWidth + (x) + 1 ) * MEM_OFFSET]

#define STOP         0
#define UP           4
#define LEFT         8
#define CROSSWISE   12
#define DIRECTION   12
#define CONTIN_UP    2
#define CONTIN_LEFT  1
#define ELEMENT     15
#define ININTI      3

//calculate_cigar_2( int size, int4 * result, char * cigar,int * cigar_int,int * direction) //,
__global__ void Calculate_ben_back(int size, short4 * result, char * cigar, int *cigar_int,  unsigned int* back, short rowWidth)
{

	int offset=blockIdx.x*blockDim.x+threadIdx.x;
	while(offset<size)
	{

		char * cigar_store;
		int * cigar_int_store;
		cigar_store=(char *) (cigar+offset*sizeof(char)*128);
		cigar_int_store=(int *) (cigar_int+offset*128);
		int segment_length;

		//startPosA == thread number within whole grid   
		int startPosA = offset;

		short4 myMaxXY = result[startPosA];
		short2 indexXY;
		indexXY.x=myMaxXY.x;
		indexXY.y=myMaxXY.y;

		int cigar_index=0;	
		if(myMaxXY.z>0)
		{
			cigar_store[cigar_index]='S';
			cigar_int_store[cigar_index]=myMaxXY.z;
			cigar_index++;
		}
		segment_length=0;


		unsigned int back8 = BACK(indexXY.x, indexXY.y);
		back8 >>= ((8 - ((indexXY.y + 1) % 8)) % 8) * 4;

		unsigned char back1 = back8 & ELEMENT; //current element of back array
		back8 >>= 4;

		unsigned char prevDirection = ININTI;// 1100 == 12 =>crosswise
		unsigned todo;
		//back 1 is the current element of back array
		while(indexXY.x>=0 && indexXY.y>=0)//(back1 & DIRECTION) //while(direction != STOP)
		{

			if( ((prevDirection & DIRECTION) == UP) && (prevDirection & CONTIN_UP) )
			{
				todo = UP;
			}
			else if( ((prevDirection & DIRECTION) == LEFT) && (prevDirection & CONTIN_LEFT) )
			{
				todo = LEFT;
			}
			else if ((back1 & DIRECTION) == UP)
			{
				todo = UP;
			}
			else if ((back1 & DIRECTION) == LEFT)
			{
				todo = LEFT;
			}
			else //if (back1 & DIRECTION == CROSSWISE)
			{
				todo = CROSSWISE;
			}

			if(prevDirection==ININTI) prevDirection=todo;
			if((prevDirection & DIRECTION)==todo) 
			{
				segment_length++;
			}
			else
			{
				//printf("             prevDirectio= %d  todo=%d\n",prevDirection,todo);
				//if(prevDirection==LEFT);
				cigar_store[cigar_index]=(prevDirection & DIRECTION);//'D';  //I D M????????
				//if(prevDirection==UP)
			//	cigar_store[cigar_index]=UP;//'I';  //I D M????????
				//if(prevDirection==CROSSWISE)
			//	cigar_store[cigar_index]=CROSSWISE;//'M';  //I D M????????

				cigar_int_store[cigar_index]=segment_length;
				cigar_index++;
				segment_length=1;
				prevDirection=todo;
			}

			if (todo == LEFT)
			{
				indexXY.x--;
				back8 = BACK(indexXY.x, indexXY.y);
				back8 >>= ((8 - ((indexXY.y + 1) % 8)) % 8) * 4; //because of the last row of back array
			}
			else if (todo == UP)
			{
				indexXY.y--;
				if((indexXY.y % 8) == 7)
					back8 = BACK(indexXY.x, indexXY.y);  //since up direction, 8 elements stored in the same int.
			}
			else //if (todo == CROSSWISE)
			{
				indexXY.x--;
				indexXY.y--;

				back8 = BACK(indexXY.x, indexXY.y);
				back8 >>= ((8 - ((indexXY.y + 1) % 8)) % 8) * 4; //because of the last row of back array
			}

			prevDirection = todo | back1&3;   //Here, back1 is used to calculate preDirection.
			//printf("prevDirection=%d  %d %d \n",prevDirection,indexXY.x,indexXY.y);
			back1 = back8 & ELEMENT;
			back8 >>= 4;
		}

		//maybe S
		//**********

		cigar_store[cigar_index]=todo;
		cigar_int_store[cigar_index]=segment_length;
		cigar_index++;
	//	printf("%d\n",cigar_index);
		if(indexXY.x>=0) 
		{	
			cigar_store[cigar_index]='S';
			cigar_int_store[cigar_index]=indexXY.x+1;
			cigar_index++;
		}	

		myMaxXY.z=indexXY.x+1;
		myMaxXY.w=cigar_index;
		result[offset]=myMaxXY;
		offset+=blockDim.x*gridDim.x;
	}

}


struct InputData
{
	char read_base[600];
	char reference_base[600];
};


int main(int artc, char* args[])
{
	uint64_t total_size=0;
	FILE * file;
        file=fopen(args[1],"r");
        int size;
	double  computation_time=0;//total_time=0;
	timespec start,finish;	 

	 char data[200][1000];   //* Here, we read in 100 sequences.
     for(int i=0;i<1;i++)
     {
        fscanf(file,"%s ", data[i]);
 }
      int row=atoi(args[2]);
      int col=atoi(args[3]);
      size=row*col;
      for(int ww=0;ww<1;ww++)
      {       
      	//Here is the sequences pairs.
      	int index=0;
        InputData * inputdata=(InputData* )malloc(100*(sizeof(InputData)));
        for(int i=0;i<1;i++)
        for(int j=0;j<1;j++)
        {
        strcpy(inputdata[index].reference_base,data[i]);
        strcpy(inputdata[index].read_base,data[j]);
  
//	printf("%s\n",inputdata[index].reference_base);
//	printf("%s\n",inputdata[index].read_base);
      index++;
        }

	 for(int j=1;j<99;j++)
        {
        strcpy(inputdata[j].reference_base,inputdata[0].reference_base);
        strcpy(inputdata[j].read_base,inputdata[0].read_base);
	}
       size=100;
		//data preparation.
		//we put all the sequence pairs into a char* array
		char * data_h_total=(char*)malloc(size * 640* sizeof (char)*2+(size*sizeof(NUM_ADD)+127)/128*128);
		NUM_ADD * data_num_add=(NUM_ADD *) (data_h_total);
		char * data_h=data_h_total+(size*sizeof(NUM_ADD)+127)/128*128;  //.thus we donot need to worry about align
		int data_size=0;
		char * data_d_total;		
		cudaMalloc( (char **) &data_d_total, (size*sizeof(NUM_ADD)+127)/128*128+size *( 640 )* sizeof (char)*2+sizeof(int)*size*4);
	//	printf("total size=%d\n",(size*sizeof(NUM_ADD)+127)/128*128+size *( 640 )* sizeof (char)*2+sizeof(int)*size*4);

		short * result_h=(short*) malloc(sizeof(short)*size*4);
	 //	printf("%d\n",sizeof(short)*size*4);
		char * cigar_h=(char *) malloc(sizeof(char)*size*128);   //Here the length of alignment is 128
		int * cigar_int_h=(int *) malloc(sizeof(int)*size*128);  //Here the length of alignment is 128

		for(int i=0;i<size;i++)
		{
			char4 reference_tep[150];
			int read_len=strlen(inputdata[i].read_base);
			int ref_len=strlen(inputdata[i].reference_base);
			int new_len=(ref_len+4-1)/4;
			total_size+=ref_len*read_len; 
			//printf("i=%d total_size=%d",i,total_size);
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
			data_h+=(read_len+128-1)/128*128;
			data_size+=(read_len+128-1)/128*128;

			memcpy(data_h,reference_tep,sizeof(char4)* new_len);
		        data_h+=(new_len*sizeof(char4)+127)/128*128;
		        data_size+=(new_len*sizeof(char4)+127)/128*128;
		}
		
		cudaError_t err;		
		int data_size_to_copy=data_size+(size*sizeof(NUM_ADD)+127)/128*128;

		cudaMemcpy(data_d_total,data_h_total,data_size_to_copy,cudaMemcpyHostToDevice);
		NUM_ADD * num_add_d=(NUM_ADD *) (data_d_total);
		char * data_d=data_d_total+(size*sizeof(NUM_ADD)+127)/128*128;
		short4 * result_d=(short4 *) (data_d_total+data_size_to_copy);
		 //printf("data size to copy=%d\n",data_size_to_copy);
		int blocksize=128;
		dim3 block(blocksize);
		dim3 grid((size+blocksize-1)/blocksize); //size/blocksize


		char * cigar;
		err=cudaMalloc( (char **) &cigar, size * (128* sizeof (char)+128*sizeof(int)));
	       if (err != cudaSuccess)
        {
printf("1   1 1 1 %s", cudaGetErrorString(err));
        }			
		int * cigar_int=(int *) (cigar+size*128*sizeof(char));
		unsigned int * direction;
		int2 * AF_maxXY;
		err=cudaMalloc((int2 **)& AF_maxXY, 640*sizeof(int2)*(size+blocksize-1)/blocksize*blocksize);// vector
		       if (err != cudaSuccess)
        {
	printf("2    23     %s", cudaGetErrorString(err));
        }


	//	cudaMalloc( (unsigned int **) & direction, size * (640*640* sizeof (unsigned int)));
	//	cudaMalloc( (unsigned int **) & direction, 640*640* sizeof (unsigned int)*(size+blocksize-1)/blocksize*blocksize);
		err=cudaMalloc( (unsigned int **) & direction, 640*(640/8)* sizeof (unsigned int)*(size+blocksize-1)/blocksize*blocksize);
		if (err != cudaSuccess)
	{
	
	printf("3      %s", cudaGetErrorString(err));
	}
		
		Calculate_ben<<<grid,block>>> (size,data_d,num_add_d,result_d,AF_maxXY, direction, 640); //result
       //Calculate_ben(int size, char *data, NUM_ADD * num_add, short4 *result, int2* AF_maxXY, unsigned int* back,  short rowWidth)
		cudaDeviceSynchronize();
       // Calculate_ben_back(int size, short4 * result, char * cigar, int *cigar_int,  unsigned int* back, short rowWidth)
		clock_gettime(CLOCK_MONOTONIC_RAW,&start);
	  Calculate_ben_back<<<grid,block>>> (size,result_d,cigar,cigar_int,direction,640); //result

	//	printf("%d\n", size*sizeof(short4));
		cudaDeviceSynchronize();
	//	cudaMemcpy(result_h,result_d,size*sizeof(short4),cudaMemcpyDeviceToHost);
	//	cudaMemcpy(cigar_h,cigar,128*sizeof(char)*size, cudaMemcpyDeviceToHost);
	//	cudaMemcpy(cigar_int_h,cigar_int,128*sizeof(int)*size,cudaMemcpyDeviceToHost);		
		
	clock_gettime(CLOCK_MONOTONIC_RAW,&finish);
        computation_time+=diff(start,finish);
	
/*		for(int i=0;i<size;i++)
		{
	//		printf("%d\n",result_h[i*4]);
	//		printf("%d\n",result_h[i*4+1]);
	//		printf("%d\n",result_h[i*4+2]);
	//		printf("%d\n",result_h[i*4+3]);
			printf("[");
			for(int j=result_h[i*4+3]-1;j>=0;j--)
			{
			printf("%d",cigar_int_h[128*i+j]);
			if(cigar_h[128*i+j]==UP)
			printf("%c",'D');
			if(cigar_h[128*i+j]==LEFT)
			printf("%c",'I');
			if(cigar_h[128*i+j]==CROSSWISE)
			printf("%c",'M');			
			if(cigar_h[128*i+j]=='S')
			printf("%c",'S');			
			if(j!=0) printf(", ");
			}
			printf("]\n");
		}
*/
        cudaFree(AF_maxXY);
		cudaFree(direction);
		free(data_h_total);
		cudaFree(data_d_total);
		free(inputdata);
		cudaFree(cigar);
		free(cigar_int_h);
		free(cigar_h);
 //     fscanf(file,"%d",&size);
        }

 //	printf(" computation_time= %e  total_time=%e \n",computation_time,0);
printf(" computation_time= %e  %d GCUPs=%lf\n",computation_time,total_size,( total_size)/computation_time/1000000000);



        return 0;
}


#undef STOP
#undef UP
#undef LEFT
#undef CROSSWISE
#undef DIRECTION
#undef CONTIN_UP
#undef CONTIN_LEFT
#undef ELEMENT



