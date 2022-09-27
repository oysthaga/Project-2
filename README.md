# Project-2

To compile Code1.cpp: 

  g++ -c Code1.cpp
  
To link: 

  g++ Code1.o -larmadillo -o Code1.exe
  
To run: 

  ./Code1.exe [Some Integer] 
  
The code requires an integer to be provided by the user. This integer is the matrix size N. 

I have run the Python scripts (Code2.py, Code3.py, Code3b.py) in the IDE Spyder by simply hitting the Run button. Code3.py assmues a 9X9 matrix, so Code1.cpp must be run first with the input 9. Similarly, Code3b.py must be run after Code1.cpp with input 99.  
