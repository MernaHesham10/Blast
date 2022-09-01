import pandas as pd

protiens = ['A','R','N','D','C','Q','E',
            'G','H','I','L','K','M','F',
            'P','S','T','W','Y','V','B','Z','X','*']

blosum62Data = pd.read_csv("BLOSUM62.csv")
blosum62Data = blosum62Data.rename(index = lambda s: protiens[int(s)])
blosum62Map = {}

for i in protiens:
    for j in protiens:
        blosum62Map[(i, j)] = blosum62Data.at[i, j]
del blosum62Data

def FindWords(inputSequence, wordLength):
    wordsList = []
    for i in range(0, len(inputSequence) - wordLength - 1):
        outputWord = inputSequence[i: i + wordLength]
        wordsList.append((outputWord, i))
    return wordsList

def FindNeighbors(words, wordLength, threshold):
    seedsList = []
    for (word, wordIndex) in words:
        wordTemp, wordLen = word, len(word)
        for i in range(wordLen):
            for aminoAcid in protiens:
                 wordTemp = wordTemp[:i] + aminoAcid + wordTemp[i + 1:]
                 finalScore = 0
                 for j in range(wordLength):
                     finalScore = finalScore + blosum62Map[(word[j], wordTemp[j])]
                 if finalScore >= threshold:
                     seedsList.append((wordTemp, word, finalScore, wordIndex))
    return seedsList

def FindHits(protien, seeds, wordLength):
    hitsList = []
    for i in range(len(seeds)):
        (wordTemp, _, _, _) = seeds[i]
        for j in range(len(protien)):
            if protien[j:j + wordLength] == wordTemp:
                hitsList.append((i, j))
    return hitsList

def HSPExtend(hits, protien, wordLength):
    hspList = []
    for i in range(len(hits)):
        seedIndex = hits[i][0]

        leftExt = hits[i][1]
        rightExt = hits[i][1] + wordLength - 1
        leftQuery = seeds[seedIndex][3]
        rightQuery = seeds[seedIndex][3] + wordLength - 1

        finalScore = seeds[seedIndex][2]
        scoreX = finalScore

        while(scoreX - finalScore < threshold):
            if(leftExt >= 0 and leftQuery >= 0):
                finalScore += blosum62Map[(protien[leftExt], inputSequence[leftQuery])]
                leftQuery -= 1
                leftExt -= 1
            if(rightExt < len(protien) and rightQuery < len(inputSequence)):
                finalScore += blosum62Map[(protien[rightExt], inputSequence[rightQuery])]
                rightExt += 1
                rightQuery += 1
            if leftExt < 0 or leftQuery < 0:
                if rightExt >= len(protien) or rightQuery >= len(inputSequence):
                    break
            if(scoreX < finalScore):
                scoreX = finalScore
        hspList.append([rightExt - 1, leftExt + 1, rightQuery - 1, leftQuery + 1, finalScore])
    return hspList

def CheckOverLapping(hsp, protien):
    NewHSPList = []
    HSP_Index = 0

    while(HSP_Index < len(hsp)):
        frExt = hsp[HSP_Index][0]
        srExt = hsp[HSP_Index + 1][0]

        flExt = hsp[HSP_Index][1]
        slExt = hsp[HSP_Index + 1][1]

        if frExt >= slExt:
            frQuery = hsp[HSP_Index][2]
            srQuery = hsp[HSP_Index + 1][2]

            flQuery = hsp[HSP_Index][3]
            slQuery = hsp[HSP_Index + 1][3]

            if frQuery >= slQuery:
                if flQuery == slQuery and frQuery == srQuery:
                    NewHSPList.append(hsp[HSP_Index])
                else:
                    if flQuery == slQuery:
                        if frQuery > srQuery:
                            NewHSPList.append(hsp[HSP_Index])
                        else:
                            NewHSPList.append(hsp[HSP_Index + 1])
                    else:
                        nSc = hsp[HSP_Index][4] + hsp[HSP_Index + 1][4]
                        counter = slExt
                        for n in range(slQuery, frQuery):
                            nSc -= blosum62Map[protien[counter]][inputSequence[n]]
                            counter += 1
                        NewHSPList.append([srExt, flExt, srQuery, flQuery, nSc])
                HSP_Index += 1
            else:
                NewHSPList.append(hsp[HSP_Index])
        else:
            NewHSPList.append(hsp[HSP_Index])
        HSP_Index += 1
    return NewHSPList

def Display(NewHSP, protien):
    for z in range(len(NewHSP)):
        x = NewHSP[z][0]
        y = NewHSP[z][1]

        print("Protien : ", protien [y : x+1])
        print("Query : ", inputSequence[NewHSP[z][3]:NewHSP[z][2]+1])
        print("Final Score : ", NewHSP[z][4])


inputSequence = 'PQGEFG'
protien = 'PYGPWGETPCGWFRRQGEHADGKRGEFQWEAAAAAPWG'

wordLen = 3
threshold = 10

words = FindWords(inputSequence, wordLen)
seeds = FindNeighbors(words, wordLen, threshold)

hits = FindHits(protien, seeds, wordLen)
HSP = HSPExtend(hits, protien, wordLen)
final = CheckOverLapping(HSP, protien)
Display(final, protien)