from Bio import Entrez
from urllib.error import HTTPError
import time
from Bio import Entrez
import nltk
import re
from nltk.tokenize import sent_tokenize
from spacy.en import English
import datetime
st=time.time()
nlp = English()
#print(time.time()  - st)
curatedPapers=[]

def getNationality(abs):
    places=[]
    for ent in abs.ents:
        if (ent.label == 347 ):
            if ((str(ent) not in places) and len(str(ent)) > 3 ):
                places.append(str(ent))

    return ' '.join(places)

#actual one used for group Size
def getcohortSzies(abs):
    cohort = ['', '']
    numModifiers = []
    ControlKeywords = 'individuals control controls normal healthy'.split()
    Patientkeywords = 'individuals patients patient case cases women men subjects participants'.split()
    for tok in abs:
        #print(tok.dep_)
        if (tok.dep_ == 'nummod'):
            #print(tok.dep, tok.orth_, tok.head.text,tok.nbor(1).text,tok.nbor(2).text)
            try:
                numModifiers.append([tok.orth_, tok.head.text + tok.nbor(1).text + tok.nbor(2).text])
            except:
                #numModifiers.append([tok.orth_, tok.head.text + tok.nbor(1).text + tok.nbor(2).text])
                pass
            ##print(numModifiers)
    #print(numModifiers)

    # controlKeywords
    for temp in numModifiers:
        temp[1] = temp[1].replace('-', '')
        if (any(word in temp[1] for word in ControlKeywords)):
            try:
                temp[0]=temp[0].replace(',', '')
                x=int(temp[0])
                cohort[1] = temp[0]
                #print("in control ", temp, cohort)
            except:
                pass

        elif (any(word in temp[1] for word in Patientkeywords)):
            try:
                temp[0] = temp[0].replace(',', '')
                x = int(temp[0])
                cohort[0] = temp[0]
                #print("in patient ", temp, cohort)
            except:
                pass
            ##print(temp)

    #print(cohort)
    return cohort


#function
def getPvaluess(abstractText):
    ###print(abstractText)
    sent_tokenize_list = sent_tokenize(abstractText)
    #sent= abstractText.sents
    PSNPList = []
    sentences=[]

    for sentence in sent_tokenize_list:
        sentences.append(sentence)
        #sentence=str(sentence)
        senten=sentence
        sentence = ''.join(c for c in sentence if c not in '(){} []')
        sentence = sentence.replace('odds ratio of','OR=').replace('odds ratio','OR').replace('-value', '').replace('value ', '').replace('meta','').replace('combined','').replace('combine','').replace('add','').replace('p of','p =')
        sentence = sentence.replace('⁻', '-').replace('¹', '1').replace('²', '2').replace('³', '3').replace('⁴', '4').replace('⁵', '5').replace('⁶', '6').replace('⁷', '7').replace('⁸', '8').replace('⁹', '9')
        #print(sentence)
        OR = re.finditer("OR[\s]?[=><≥≤]?[ \s]?[-+]?\d[\.·]?\d*[\s]*", sentence)
        if re.findall("[pP][\s]?[=><≥≤][ \s]?[-+]?\d[\.·]?\d*[\s]?[x×]?[\s]?[0-9eE]?[0-9]?[\s]?[+-⁻]?[\s]?..", sentence,flags=re.IGNORECASE):#resultsSentence

            # print(sentence)
            sentence = re.sub('(?=)(\s)(?=[\d+=><≥≤×*x(10)-])', '', sentence)
            sentence = re.sub('([×*xX]10)', 'E', sentence)
            sentence = sentence.replace('=10', '=1e').replace('>10','>1e' ).replace('<10','<1e' ).replace('≥10','≥1e').replace('≤10','≤1e')
            sentence = re.sub('(\s?[×*xX]\s?E)', 'E', sentence)
            #print(sentence)
            SNPList   = []
            PvalList = []
            ORList=[]
            #P = re.finditer("[pP][\s]?[=><≥≤][ \s]?[-+]?\d[\.]?\d*[\s]?[x×]?[\s]?([eE][+-]\d*)?", sentence)

            P = re.finditer("[pP][\s]?[=><≥≤][ \s]?[-+]?(\d[\.·]?\d*)?[\s]?[x×]?[\s]?([eE][⁻+-]\d*)?[¹²³⁴⁵⁶⁹]*", sentence)

            for m in P:
                temp=''.join(c for c in m.group(0) if c not in 'pP><≤≥= ')
                PvalList.append([ temp,m.start()])

            #print(PvalList)
            for m in OR:
                temp = ''.join(c for c in m.group(0) if c not in 'OR><≤≥= ')
                ORList.append([temp, m.start()])
            #L=re.findall("[pP][\s]?[=><≥≤][ \s]?[-+]?\d[\.]?\d*[\s]?[x×]?[\s]?[0-9eE]?[0-9]?[\s]?[+-⁻]?[\s]?..", sentence,flags=re.IGNORECASE)#resultsSentence)
            S = re.finditer("(rs\d[-]*[\w]*)",sentence)

            for m in S:
                SNPList.append([m.group(), m.start()])

            #print(sentence)
            #print(SNPList,PvalList)
            #ignore pValues before SNP mentions_RARE
            if (len(PvalList) > len(SNPList) and SNPList ):
                if (len(SNPList)==1 ):
                    for temp in range (1,len(PvalList)):
                        SNPList.append(SNPList[0])


                else:
                    while (PvalList and (SNPList[0][1]-PvalList[0][1] >0)):
                        PvalList.pop(0)

            elif (len (PvalList) ==1 and len(SNPList)>1):
                for temp in range(1, len(SNPList)):
                    PvalList.append(PvalList[0])

            if (len(ORList) > len(SNPList) and SNPList):
                while (ORList and (SNPList[0][1] - ORList[0][1] > 0)):
                    ORList.pop(0)

            for i in SNPList:
                x=['-1','']
                y=['-1','']

                if PvalList:
                    x=PvalList.pop(0)
                if ORList:
                    y=ORList.pop(0)

                if (x[0]!='-1' and y[0]!='-1'):
                    PSNPList.append([i[0], x[0],y[0], senten])
                elif (x[0]=='-1' and y[0]!='-1'):
                    PSNPList.append([i[0], '', y[0], senten])
                elif (x[0] != '-1' and y[0] == '-1'):
                    PSNPList.append([i[0], x[0], '', senten])

                ##print([i[0], x[0],y[0], sentence])
    return PSNPList


def getPvaluesss(abstractText):
    ###print(abstractText)
    sent_tokenize_list = sent_tokenize(abstractText)
    #sent= abstractText.sents
    PSNPList = []
    sentences=[]
    for sentence in sent_tokenize_list:
        sentences.append(sentence)
        #sentence=str(sentence)
        senten=sentence
        sentence = ''.join(c for c in sentence if c not in '(){} []')
        sentence = sentence.replace('-value', '').replace('value ', '').replace('meta','').replace('combined','').replace('combine','').replace('add','')
        sentence = sentence.replace('⁻', '-').replace('¹', '1').replace('²', '2').replace('³', '3').replace('⁴', '4').replace('⁵', '5').replace('⁶', '6').replace('⁷', '7').replace('⁸', '8').replace('⁹', '9')
        #print(sentence)
        if re.findall("[pP][\s]?[=><≥≤][ \s]?[-+]?\d[\.·]?\d*[\s]?[x×]?[\s]?[0-9eE]?[0-9]?[\s]?[+-⁻]?[\s]?..", sentence,flags=re.IGNORECASE):#resultsSentence
            #print(sentence)
            sentence = re.sub('(?=)(\s)(?=[\d+=><≥≤×*x(10)-])', '', sentence)
            sentence = re.sub('([×*xX]10)', 'E', sentence)
            sentence = sentence.replace('=10', '=1e').replace('>10','>1e' ).replace('<10','<1e' ).replace('≥10','≥1e').replace('≤10','≤1e')
            sentence = re.sub('(\s?[×*xX]\s?E)', 'E', sentence)
            #print(sentence)
            SNPList = []
            PvalList = []
            ORList=[]
            #P = re.finditer("[pP][\s]?[=><≥≤][ \s]?[-+]?\d[\.]?\d*[\s]?[x×]?[\s]?([eE][+-]\d*)?", sentence)

            P = re.finditer("[pP][\s]?[=><≥≤][ \s]?[-+]?(\d[\.·]?\d*)?[\s]?[x×]?[\s]?([eE][⁻+-]\d*)?[¹²³⁴⁵⁶⁹]*", sentence)
            OR = re.finditer("OR[\s]?[=><≥≤][ \s]?[-+]?\d[\.·]?\d*", sentence)
            for m in P:
                temp=''.join(c for c in m.group(0) if c not in 'pP><≤≥= ')
                PvalList.append([ temp,m.start()])

            for m in OR:
                temp = ''.join(c for c in m.group(0) if c not in 'OR><≤≥= ')
                ORList.append([temp, m.start()])
            #L=re.findall("[pP][\s]?[=><≥≤][ \s]?[-+]?\d[\.]?\d*[\s]?[x×]?[\s]?[0-9eE]?[0-9]?[\s]?[+-⁻]?[\s]?..", sentence,flags=re.IGNORECASE)#resultsSentence)
            S = re.finditer("(rs\d[-]*[\w]*)",sentence)

            for m in S:
                SNPList.append([m.group(), m.start()])
            #print(sentence)
            #print(SNPList,PvalList)
            #ignore pValues before SNP mentions_RARE
            if (len(PvalList) > len(SNPList) and SNPList ):
                while (PvalList and (SNPList[0][1]-PvalList[0][1] >0)):
                    PvalList.pop(0)

            if (len(ORList) > len(SNPList) and SNPList):
                while (ORList and (SNPList[0][1] - ORList[0][1] > 0)):
                    ORList.pop(0)

            for i in SNPList:
                x=['','']
                y=['','']

                if PvalList:
                    x=PvalList.pop(0)
                if ORList:
                    y=ORList.pop(0)
                PSNPList.append([i[0], x[0],y[0], senten])
                ##print([i[0], x[0],y[0], sentence])
    return PSNPList

def getPubMedId(targetD,year1,year2,nPub):
    #targetD = "\""+targetD + "\" AND (\"genome-wide\" OR \"genome AND identification\" OR \"genome AND association\" OR \"candidate AND gene\" OR \"SNP\") "
    #targetD = targetD + "AND Genetic Predisposition to Disease"
    #targetD=targetD + " AND ( genome-wide OR (genome AND identification) OR (genome AND association) OR (candidate AND gene) OR (SNP) ) "
    #targetD = targetD + " AND (Polymorphism, Single Nucleotide OR genetic Predisposition to Disease)"

    #targetD= "\""+ targetD+"\"[title] AND ((\"polymorphism, single nucleotide\"[MeSH Terms]) OR(\"genetic predisposition to disease\"[MeSH Terms]))"
    #targetD="("+targetD+") AND ((\"polymorphism, single nucleotide\"[MeSH Terms]) OR(\"genetic predisposition to disease\"[MeSH Terms]))"
    #targetD = "\"" + targetD + "\" AND ((\"polymorphism, single nucleotide\"[MeSH Terms]) OR(\"genetic predisposition to disease\"[MeSH Terms]))"

    #targetD= "\"" + targetD + "\" AND ( (\"polymorphism, single nucleotide\"[MeSH Terms]) OR(\"genetic predisposition to disease\"[MeSH Terms]) OR \"genome-wide\" OR \"genome AND identification\" OR \"genome AND association\" OR \"candidate AND gene\" OR \"SNP\")"
    #targetD = targetD + "[All Fields] AND(\"Polymorphism, Single Nucleotide\"[Mesh Terms] OR \"Genetic Predisposition to Disease\"[Mesh Terms] OR \"Genome-Wide Association Study\"[Mesh Terms] OR \"genome-wide\" OR \"genome AND identification\" OR \"genome AND association\" OR \"candidate AND gene\")"

    targetD ="(("+"\""+year1+"\"[Date - Publication] : "+"\""+year2+"\"[Date - Publication]))"+"AND "+"\"" + targetD + "\" " +"AND(\"Polymorphism, Single Nucleotide\"[Mesh Terms] OR \"Genetic Predisposition to Disease\"[Mesh Terms] OR \"Genome-Wide Association Study\"[Mesh Terms])"
    Entrez.email = "n.s.tawfik@uu.nl"  # Always tell NCBI who you are
    print("Searchig PubMed IDs relevant to keyword: ",targetD)
    handle = Entrez.esearch(db="pubmed", term=targetD,RetMax=nPub)
    Rlist = Entrez.read(handle)
    print("The total number of Extracted Publications is " + Rlist["Count"])
    return Rlist["IdList"]

def getPubMedArticle(IDs):
    #print(len(IDs),IDs)
    AllRecords = []
    #divide Idlist into smaller lists to avoid 403 **BOT Issues**
    IDlists=[IDs[i:i + 200] for i in range(0, len(IDs), 200)]
    for IDlist in IDlists:
        attempt=1
        while attempt <= 3:
            try:
                fetch_handle = Entrez.efetch(db="pubmed", id=IDlist, rettype="null", retmode="xml", parsed=True)
                if 'fetch_handle' in locals():
                    data = Entrez.read(fetch_handle)
                    fetch_handle.close()
                    for i in data['PubmedArticle']:
                        AllRecords.append([i])
                    #print(i)
                attempt = 4
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print("Received error from server %s" % err)
                    print("Attempt %i of 3" % attempt)
                    attempt += 1
                    time.sleep(15)
                else:
                    raise
    #print(len(AllRecords))
    return AllRecords


def extractInfo(PaperRecords):
    curatedPapers=[]
    j=0
    numberOfAbsWithSNP = 0
    for record in PaperRecords:
        # getabstract
        try:  # check for abstract
            abstract = record[0]['MedlineCitation']['Article']['Abstract']['AbstractText']
        except KeyError as e:
            continue
        abstractText = ""
        for k in abstract:
            abstractText += (str(k))
        sent = nltk.tokenize.wordpunct_tokenize(abstractText)
        SNPstr = ([w for w in set(sent) if w.startswith('rs')] or [w for w in set(sent) if w.startswith('ss')])
        numberOfAbs=len(PaperRecords)
        if SNPstr:
            numberOfAbsWithSNP+=1
            # spaCy
            abstractText = abstractText.replace('(', '').replace(')', '').replace(' n ', 'with n').replace('p.', 'p')
            abstractTextt = nlp(abstractText)
            #print(record[0]['MedlineCitation']['PMID'])

            pValuesExtracted = getPvaluess(abstractText)
            if(pValuesExtracted):
                #print(pValuesExtracted)
                population = getNationality(abstractTextt)
                gSize = getcohortSzies(abstractTextt)
            # gSize=getGroupSize(abstractText)
                patientSize = gSize[0]
                controlSize = gSize[1]
                #print(patientSize, controlSize)
                writeInfo(j, record, abstractText, population, patientSize, controlSize, SNPstr, pValuesExtracted)
                j += 1

    return numberOfAbs, numberOfAbsWithSNP


def writeInfo(j,record,abstractText,population,patientSize,controlSize,SNPstr, pValuesExtracted):
    curatedPapers.append([])
    curatedPapers[j].append(str(record[0]['MedlineCitation']['PMID']))# ID
    curatedPapers[j].append(str(record[0]['MedlineCitation']['Article']['ArticleTitle']))  # title
    curatedPapers[j].append(record[0]['MedlineCitation']['Article']['ArticleDate']) # date
    curatedPapers[j].append(str(abstractText))  # abstract
    curatedPapers[j].append(population)  # population
    curatedPapers[j].append(patientSize)  # Patient Group Size
    curatedPapers[j].append(controlSize)  # Control Group Size
    curatedPapers[j].append(str(SNPstr))  # SNP list found
    curatedPapers[j].append(pValuesExtracted) #extracted P values
    curatedPapers[j].append(0)  #extracted frequency

    if curatedPapers[j][2]:
        x = curatedPapers[j][2][0].get('Day') + curatedPapers[j][2][0].get('Month') + curatedPapers[j][2][0].get('Year')
        x = curatedPapers[j][2][0].get('Day').zfill(2) + curatedPapers[j][2][0].get('Month').zfill(2) + curatedPapers[j][2][0].get('Year')
        try:
            curatedPapers[j][2] = datetime.datetime.strptime(x, "%d%b%Y").date()
        except:
            try:
                curatedPapers[j][2] = datetime.datetime.strptime(x, "%d%m%Y").date()
            except:
                print(curatedPapers[j][2],"ERROR")

    else:
        curatedPapers[j][2] =''


def curate(disease,year1,year2,nPub):
    f = open('myfile.txt', 'w')
    f.write(disease)
    f.flush()
    curatedPapers.clear()
    if '-' in disease:
        disease=disease +" OR "+ disease.replace("-", "")
    st = time.time()
    IDlist=getPubMedId(disease,year1,year2,nPub)
    f.write("\nTime to get IDs\n",)
    #print(len(IDlist),IDlist)
    f.write(str(time.time() - st))
    f.flush()
    st = time.time()
    PaperRecords=getPubMedArticle(IDlist)
    f.write("\nTime to download records")
    f.write(str(time.time() - st))
    f.flush()
    st = time.time()
    nAbs,nSP=extractInfo(PaperRecords)
    #getFrequency()
    #print(nAbs,nSP)
    f.write("\nInformation extraction time")
    f.write(str(time.time() - st))
    f.flush()
    print(nAbs, nSP)
    return(curatedPapers,nAbs,nSP)

#curate('SKIN PIGMENTATION','2000','2017', 3000)