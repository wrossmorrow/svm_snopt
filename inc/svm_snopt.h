/*
 *  svm_snopt.h
 *  
 *
 *  Created by W. Ross Morrow on 10/23/13.
 *  Copyright 2013 Iowa State University. All rights reserved.
 *
 */

#ifndef SVMSNOPT_H
#define SVMSNOPT_H

int svm_snopt(int   extA ,      // number of attributes
              int * extL ,      // A-element vector of number of levels per attribute
              int   extJ ,      // number of profiles
              int * extX ,      // J x A element matrix of profiles (level indexing form)
              int * exty ,      // J-vector of classification responses (-1,1)
              int   extC ,      // soft-margin penalty term
              char  ic ,        // initial condition code
              double * x0 ,     // optional prescribed initial condition for parameters
              char prob ,       // problem code ('d' for "dual", anything else for primal)
              double * w ,      // computed parameters: actual parameters
              double * b ,      // computed parameters: constant term
              char * optfn ,    // options filename
              char * logfn);    // log file name

#endif