/*
 *  csvreader.h
 *  
 *
 *  Created by W. Ross Morrow on 7/2/13.
 *  Copyright 2013 Iowa State University. All rights reserved.
 *
 */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * GENERIC CSV Reader for coefficient data * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void CSVReader(char * fnm , 
			   int M , 
			   int N , 
			   double * data , 
			   enum CBLAS_ORDER order );