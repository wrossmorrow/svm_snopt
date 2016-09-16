/*
 *  csvreader.c
 *  
 *
 *  Created by W. Ross Morrow on 7/2/13.
 *  Copyright 2013 Iowa State University. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vecLib/cblas.h>

#include "csvreader.h"

#define BUFFER_SIZE 1024

void CSVReader(char * fnm , 
			   int M , 
			   int N , 
			   double * data , 
			   enum CBLAS_ORDER order )
{
	
	char buf[BUFFER_SIZE];
	int i = 0, j;
	FILE * Data;
	char * tok;
	
	if( fnm == NULL ) {
		printf("CSV Reader Error - No valid filename given.\n");
		return;
	}
	
	if( M <= 0 || N <= 0 ) {
		printf("CSV Reader Error - Given sizes (%i,%i) not valid; both must be positive.\n",M,N);
		return;
	}
	
	if( data == NULL ) { 
		printf("CSV Reader Error - Array not provided. Cannot store data.\n");
		return;
	}
	
	Data = fopen(fnm,"r");
	if( Data == NULL ) { printf("CSV Reader Error - cannot read file %c%s%c.\n",'"',fnm,'"'); }
	else {
		
		switch( order ) {
				
			case CblasColMajor: // column-major storage
				
				i = 0;
				
				while( fgets(buf, sizeof(buf), Data) != NULL && i < M )	{
					
					j = 0;
					tok = strtok( buf, "," );
					
					while( tok != NULL && j < N ) {
						
						// read in coefficient data, storing
						// coefficients as the - rows - of a 
						// M X N matrix in row-major order
						data[ M * j + i ] = atof( tok );
						
						// get next CSV entry
						tok	= strtok( NULL , "," );
						
						// increment column counter
						j++;
						
					}
					
					// Take care of case where no more tokens remain, but we haven't
					// read enough entries in this row. This is somewhat extraneous, 
					// as we could just not read any further into this row of the array. 
					if( tok == NULL && j < N ) {
						
						printf("WARNING -- expected %i coefficients per row, but\n",N);
						printf("           Row %i has only %i coefficients.\n",i+1,j+1);
						printf("           Assigning remaining coefficients as zero.\n");
						
						while( j < N ) { data[ M * j + i ] = 0.0; j++; }
						
					}
					
					// increment row counter
					i++;
					
				}
				break;
				
			default: // row-major storage
				
				i = 0;
				
				while( fgets(buf, sizeof(buf), Data) != NULL && i < M )	{
					
					j = 0;
					tok = strtok( buf, "," );
					
					while( tok != NULL && j < N ) {
						
						// read in coefficient data, storing
						// coefficients as the - rows - of a 
						// M X N matrix in row-major order
						data[ N * i + j ] = atof( tok );
						
						// get next CSV entry
						tok	= strtok( NULL , "," );
						
						// increment column counter
						j++;
						
					}
					
					// Take care of case where no more tokens remain, but we haven't
					// read enough entries in this row. This is somewhat extraneous, 
					// as we could just not read any further into this row of the array. 
					if( tok == NULL && j < N ) {
						
						printf("WARNING -- expected %i coefficients per row, but\n",N);
						printf("           Row %i has only %i coefficients.\n",i+1,j+1);
						printf("           Assigning remaining coefficients as zero.\n");
						
						while( j < N ) { data[ N * i + j ] = 0.0; j++; }
						
					}
					
					i++;
					
				}
				break;
				
		}
		
		fclose(Data);
		
	}
	
}