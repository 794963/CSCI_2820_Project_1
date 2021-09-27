  string caseName = "";
  for (int count =  0 ; count < argc; ++count)
    {
      if ( !strcmp(argv[count],"-case"  ) ) caseName  =argv[count+1];
    }
  
  if ( caseName == "" ) { cout << "Choose a code with -case <ID>\n"; exit(0); }
  
  // --------------------------
  // Initialize the array
  // --------------------------

  // 4x6
  
  int m1 = 4;
  int n1 = 6;
  double A1[m1][n1] = {  {   4. , -1. , -1. , -1. , -3. , -5.  } ,
			 { -10. ,  4. , -1. , -1. , -3. , -10. } ,
			 { -1.  , -1. ,  4. , -1. , -3. , -7.  } ,
			 { -1.  , -1. , -1. ,  4. , -6. , -2.  }    };
  double b1[m1]    =     {  1.  ,  3. ,  4. ,  2.              };
  
  // 4x4 with unique solution
  
  int   m2 = 4;
  int   n2 = 4;
  double A2[m2][n2] = {  {4.,-1.,-1.,3.} , {-1.,4.,-3.,10.} , {-1.,-1.,4.,5.} , {-3.,-1.,-1.,8. } };
  double b2[m2]    = { 1., 3., 4.,4.};
  
  // From notes (1.13a) and book, page 14
  
  int   m3 = 4;
  int   n3 = 4;
  double A3[m3][n3] = {   {  0. , -3. , -6. ,  4. }, { -1. , -2. , -1. ,  3. } , { -2. , -3. ,  0. ,  3. }  , {  1. ,  4. ,  5. , -9. }   };
  double b3[m3]    =  { 9.  ,  1. , -1. , -7. };
  
  // From book, page 44
  
  int   m4 = 3;
  int   n4 = 3;
  double A4[m4][n4] = {   {  3. , 5. , -4. }, { -3. , -2. , 4. } , { 6. , 1. , -8.  } };
  double b4[m4]    =  { 0., 0., 0.};
  
  // From book, page 46
  
  int   m5 = 3;
  int   n5 = 3;
  double A5[m5][n5] = {   {  3. ,  5. , -4.  } ,
			  { -3. , -2. ,  4.  } ,
			  {  6. ,  1. , -8.  } };
  double b5[m5]     =     {  7. , -1. , -4.  };

  
  // 4x4 with non-unique solution
  
  int   m6 = 4;
  int   n6 = 4;
  double A6[m6][n6] = {  {4.,-1.,-1.,3.} , {3.,-10.,-1.,3.} , {-1.,-1.,-2.,-6.} , {-3.,-1.,-1.,-4.} };
  double b6[m6]     = { 1., 3., 4.,4.};

  // Exercise 12, page 48
  
  int   m7 = 4;
  int   n7 = 6;
  double A7[m7][n7] = {  { 1  ,  5  ,  2  ,  -6  ,  9  ,  0 }, 
  	               { 0  ,  0  ,  1  ,  -7  ,  4  , -8 },
		       { 0  ,  0  ,  0  ,   0  ,  0  ,  1 },
		       { 0  ,  0  ,  0  ,   0  ,  0  ,  0 }};
   double b7[m7]     = {  0  ,  0  ,   0 , 0};

  // HW (homeworks)
  
  int   m8 = 3;
  int   n8 = 4;
  double A8[m8][n8] = {  {   1  ,  -7  ,  0  ,  6 },
  	               {   0  ,  0   ,  1  , -2 },
		       {  -1  ,  7   , -4  ,  2 } };
  double b8[m8]     = {  5   ,   -3,   7  };

  // --------------------------
  // Set the case
  // --------------------------

  m = 0;
  n = 0;

  if ( caseName == "4x6"       )  {      m = m1;      n = n1;    }
  if ( caseName == "unique"    )  {      m = m2;      n = n2;    }
  if ( caseName == "page14"    )  {      m = m3;      n = n3;    }
  if ( caseName == "page44"    )  {      m = m4;      n = n4;    }
  if ( caseName == "page46"    )  {      m = m5;      n = n5;    }
  if ( caseName == "nonunique" )  {      m = m6;      n = n6;    }
  if ( caseName == "ex12"      )  {      m = m7;      n = n7;    }
  if ( caseName == "hw"        )  {      m = m8;      n = n8;    }

  if ( m == 0 ) { cout << "Bad case selection\n"; exit(0);}
  

  double b[m];
  double A[m][n];

#define iLOOP   for ( int i = 0 ; i < m ; ++i )
#define jLOOP   for ( int j = 0 ; j < n ; ++j )

  if ( caseName == "4x6"       )   iLOOP  { b[i] = b1[i];  jLOOP A[i][j] = A1[i][j];   }
  if ( caseName == "unique"    )   iLOOP  { b[i] = b2[i];  jLOOP A[i][j] = A2[i][j];   }  
  if ( caseName == "page14"    )   iLOOP  { b[i] = b3[i];  jLOOP A[i][j] = A3[i][j];   }
  if ( caseName == "page44"    )   iLOOP  { b[i] = b4[i];  jLOOP A[i][j] = A4[i][j];   }
  if ( caseName == "page46"    )   iLOOP  { b[i] = b5[i];  jLOOP A[i][j] = A5[i][j];   }
  if ( caseName == "nonunique" )   iLOOP  { b[i] = b6[i];  jLOOP A[i][j] = A6[i][j];   }
  if ( caseName == "ex12"      )   iLOOP  { b[i] = b7[i];  jLOOP A[i][j] = A7[i][j];   }
  if ( caseName == "hw"        )   iLOOP  { b[i] = b8[i];  jLOOP A[i][j] = A8[i][j];   }

   cout << "\n";
   cout << "\n";
   cout << "\n";
   cout << "========================================================\n";
   cout << "========================================================\n";
   cout << "Running case " + caseName << " \n";
   cout << "========================================================\n";
   cout << "========================================================\n";
   cout << "\n";
   cout << "\n";
   cout << "\n";


