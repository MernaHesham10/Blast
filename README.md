# Blast

GENERAL DESCRIPTION : Implementation of Blast technique to search for a protein query in a protein database.

# INPUT :

    1- A protein sequence query

    2- Word threshold

    3- Word length

    4- HSP threshold 

# OUTPUT

  1- HSP(s) with their score if they exist. Along with the ID of the sequence in the database that the 
  hsp was found in (For example, the first sequence would have ID of zero, the second sequence 
  would have ID of 1…etc)

  2- If not, then output a message indicating that the query wasn’t found. 
  SPECIFIC NOTES 

      1- You can use any similarity matrix of your choice (Ex. BLOSUM80, BLOSUM62…etc)

      2- You should have a database in the form of a text file that is being read every time the program 
      runs that contains multiple protein sequences.

  3- You should follow the steps that we took in the lab and any specifics not mentioned in the steps 
  or this document can be assumed by your team.

  4- When generating the neighboring words, you can assume you would only change one amino acid 
  at a time only. So for example, if the word is QAT, you should change Q to all its possible values 
  with AT staying the same, then change A to all its possible values with QT staying the same and 
  finally change T to all its possible values with QA staying the same

  5- You should stop extending to find the hsp when: 

    1- The query is done from both sides

    2- The score of the current alignment dropped by more than “X” from the best score that you 
    got, X is a number that you can specify. 

    For example, if X = 3, and extending your sequence had a best score of 15, then after 
    extending it more it had the score of 13, you would still continue extending till it reaches 11 
    because the difference between the current score and the best score is more than 3 (15-11 = 4)
  
  6 - If after extending, the hsp’s score is less than the HSP threshold the user specified, then ignore that alignment.
