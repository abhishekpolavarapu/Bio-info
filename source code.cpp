#include <iostream>
#include<fstream>
#include<sstream>
#include <vector>
#include<stdlib.h>
#include<iomanip>
#include<string>

using namespace std;
void extractData(string s, ofstream & outfile);
void extracthelixData(string s, ofstream & outfile);
void extractsidechain(string s, ofstream & outfile);
void extractsheetData(string s, ofstream & outfile);
void extractBetastrand(string s, ofstream & outfile);
const int MAX_VALUES = 20000;
struct atomData
{
    string atomName;
    int atomNumber;
    string atomType;
    string aminoAcidType;
    string chain;
    int residueNumber;
    double xCoordinate;
    double yCoordinate;
    double zCoordinate;
    string unknown;
    string unknown2;

};
struct helixData
{
    string recordType;
    int serialnumber;
    string identifier;
    string residueName;
    string chainIdentifier;
    int residueSequenceNumber;
    string unknown1;
    string terminalResidueName;
    string chainidentifier;
    int residueSequenceNumber2;
    string unknown2;

};


struct sheetData
{
    string recordType;
    string residueName;
    int sequenceNumber;
    string aminoType1;
    string terminalResidueName;
    string aminoType2;
    int terminalISequenceNumber;
};


std::vector<atomData> atoms;
std::vector<helixData> helix;
std::vector<sheetData> sheets;

int main()
{
    string s,filename;
    ifstream inFile,in,inhelix,insidechain;
    ofstream outFile,helixout,sidechainout,sheetout,out;


    cout<<"Enter the file name:";
    cin>>filename;
    inFile.open(filename.c_str());
    outFile.open("tempoutput.pdb");
    helixout.open("helix.temp");
    sidechainout.open("tempsidechain.pdb");
    sheetout.open("sheet.temp");
    out.open("beta.pdb");

   if(!inFile)
   {
       cout<<" Unable to open the file.Please try again!"<<endl;
       exit(0);

   }

   else
    {
        while(!inFile.eof())
    {
        getline(inFile,s);
        extractsidechain(s,sidechainout);
        extractData(s, outFile);

        extracthelixData(s, helixout);
        extractsheetData(s,sheetout);
        extractBetastrand(s,out);

    }

    }
cout<<"Processing....."<<endl;
cout<<" "<<endl;
cout<<" The Reduced representation is written to the following files : "<<endl;
cout<<"  --- sideChain-output.pdb"<<endl;
cout<<"  --- helix-output.pdb"<<endl;
cout<<"  --- betaStrand-output.pdb"<<endl;
inFile.close();
outFile.close();
helixout.close();
sidechainout.close();
in.open("tempoutput.pdb");
inhelix.open("helix.temp");
outFile.open("helix-output.pdb");
string one[MAX_VALUES],three[MAX_VALUES],four[MAX_VALUES];
char five[MAX_VALUES],last[MAX_VALUES];
double  six[MAX_VALUES]={0}, xcor[MAX_VALUES]={0}, ycor[MAX_VALUES]={0}, zcor[MAX_VALUES]={0},un1[MAX_VALUES]={0},un2[MAX_VALUES]={0},xtemp=0,ytemp=0,ztemp=0;
int index = 0,length=0 ,index1=0,length1=0,index2=0,length2=0,two[MAX_VALUES]={0};

    if(!in){
        cout<<"\n\n  Error Opening the file \n\n";
    }

    while(in>>one[index]>>two[index]>>three[index]>>four[index]>>five[index]>>six[index]>>xcor[index]>>ycor[index]>>zcor[index]>>un1[index]>>un2[index]>>last[index])
    {

       index++;
    }
    length = index;

string hone[MAX_VALUES],htwo[MAX_VALUES],hthree[MAX_VALUES],hfour[MAX_VALUES],hfive[MAX_VALUES],hseven[MAX_VALUES],height[MAX_VALUES],hten[MAX_VALUES],heleven[MAX_VALUES],htwelve[MAX_VALUES],hthreeteen[MAX_VALUES];
int hsix[MAX_VALUES]={0}, hnine[MAX_VALUES] ={0};

while(inhelix>>hone[index1]>>htwo[index1]>>hthree[index1]>>hfour[index1]>>hfive[index1]>>hsix[index1]>>hseven[index1]>>height[index1]>>hnine[index1]>>hten[index1]>>heleven[index1]>>htwelve[index1]>>hthreeteen[index1])
    {
        index1++;

    }
    length1 = index1;

    for(int j=0;j<=length1;j++)
    {
    for(int i=0;i<length-3;i++)
    {
        if((six[i]>= hsix[j] && six[i]<=hnine[j]))
       {xtemp= (xcor[i]+xcor[i+1]+xcor[i+2]+xcor[i+3])/4;
        ytemp= (ycor[i]+ycor[i+1]+ycor[i+2]+ycor[i+3])/4;
        ztemp= (zcor[i]+zcor[i+1]+zcor[i+2]+zcor[i+3])/4;
        outFile<<one[i]<<setw(7)<<right<<two[i]<<setw(3)<<right<<'S'<<setw(5)<<four[i]<<setw(2)<<'A'<<setw(4)<<six[i]<<setw(12)<<xtemp<<setw(8)<<ytemp<<setw(8)<<ztemp<<endl;

    }

    }
 outFile<<"TER"<<endl;
    }
    in.close();
    inFile.close();
    outFile.close();
////////////////////////////////////////////////////////////////////////////////////////
in.open("tempsidechain.pdb");
outFile.open("sideChain-output.pdb");
string scone[MAX_VALUES],scthree[MAX_VALUES],scfour[MAX_VALUES];
char scfive[MAX_VALUES],sclast[MAX_VALUES];
double  scsix[MAX_VALUES]={0}, scxcor[MAX_VALUES]={0}, scycor[MAX_VALUES]={0}, sczcor[MAX_VALUES]={0},scun1[MAX_VALUES]={0},scun2[MAX_VALUES]={0},xtemp1=0,ytemp1=0,ztemp1=0;
int sctwo[MAX_VALUES]={0}, count =1;


if(!in){
        cout<<"\n\n  Error opening the side chain file\n\n";
    }

 while(in>>scone[index2]>>sctwo[index2]>>scthree[index2]>>scfour[index2]>>scfive[index2]>>scsix[index2]>>scxcor[index2]>>scycor[index2]>>sczcor[index2]>>scun1[index2]>>scun2[index2]>>sclast[index2])
    {

       index2++;
    }
    length2 = index2;

    for(int k=0;k<length2;k++)
    {
        if(scsix[k]==scsix[k+1])
        {
            count++;
         xtemp1 = (xtemp1 + scxcor[k]);
         ytemp1 = (ytemp1 + scycor[k]);
         ztemp1 = (ztemp1 + sczcor[k]);
            continue;
        }

          xtemp1 = (xtemp1+scxcor[k])/count;
          ytemp1= (ytemp1+scycor[k])/count;
          ztemp1 = (ztemp1+sczcor[k])/count;
        outFile<<scone[k]<<setw(7)<<right<<sctwo[k]<<setw(3)<<right<<'S'<<setw(5)<<scfour[k]<<setw(2)<<'A'<<setw(4)<<scsix[k]<<setw(12)<<xtemp1<<setw(8)<<ytemp1<<setw(8)<<ztemp1<<endl;
        count=1; xtemp1=0; ytemp1=0; ztemp1=0;
    }
in.close();
outFile.close();
///////////////////////////////////////////////////

in.open("sheet.temp");
inFile.open("beta.pdb");
outFile.open("betaStrand-output.pdb");

string sheetone[MAX_VALUES],sheettwo[MAX_VALUES],sheetthree[MAX_VALUES],sheetfive[MAX_VALUES],sheetsix[MAX_VALUES];
int sheetfour[MAX_VALUES]={0},sheetseven[MAX_VALUES]={0},index3=0,index4=0,length3=0,length4=0;
if(!in){
        cout<<"\n\n  Error opening the Beta Strand file\n\n";
    }
while(in>>sheetone[index3]>>sheettwo[index3]>>sheetthree[index3]>>sheetfour[index3]>>sheetfive[index3]>>sheetsix[index3]>>sheetseven[index3])
{
    index3++;
}
length3 = index3;

string bsone[MAX_VALUES],bsthree[MAX_VALUES],bsfour[MAX_VALUES];
char bsfive[MAX_VALUES],bslast[MAX_VALUES];
double   bsxcor[MAX_VALUES]={0}, bsycor[MAX_VALUES]={0}, bszcor[MAX_VALUES]={0},bsun1[MAX_VALUES]={0},bsun2[MAX_VALUES]={0},xtemp2=0,ytemp2=0,ztemp2=0;
int bstwo[MAX_VALUES]={0},bssix[MAX_VALUES]={0};

while(inFile>>bsone[index4]>>bstwo[index4]>>bsthree[index4]>>bsfour[index4]>>bsfive[index4]>>bssix[index4]>>bsxcor[index4]>>bsycor[index4]>>bszcor[index4]>>bsun1[index4]>>bsun2[index4]>>bslast[index4])

{

    index4++;
}
length4 = index4;

 for(int l=0;l<=length3;l++)
    {
    for(int m=0;m<length4-3;m++)
    {
        if((bssix[m]>= sheetfour[l]) && (bssix[m]<=sheetseven[l]))
       {xtemp2= (bsxcor[m]+bsxcor[m+1]+bsxcor[m+2]+bsxcor[m+3])/4;
        ytemp2= (bsycor[m]+bsycor[m+1]+bsycor[m+2]+bsycor[m+3])/4;
        ztemp2= (bszcor[m]+bszcor[m+1]+bszcor[m+2]+bszcor[m+3])/4;
        outFile<<bsone[m]<<setw(7)<<right<<bstwo[m]<<setw(3)<<right<<'S'<<setw(5)<<bsfour[m]<<setw(2)<<'A'<<setw(4)<<bssix[m]<<setw(12)<<xtemp2<<setw(8)<<ytemp2<<setw(8)<<ztemp2<<endl;
    }
    }

    }


    return 0;
}


void extractData(string s, ofstream & outfile)
{
    string word;
    atomData aData;
    int value;
     word=s.substr(0,5);
     string temp="";
         word=s.substr(0,4);
         string word1="";
         string word2="";
        if(word=="ATOM")
        {
            word1=s.substr(13,2);
        }

        if(word=="ATOM" && word1=="CA")
        {
            outfile<<endl;
            aData.atomName=word;
            outfile<<word<<"   ";
            word=s.substr(7,4);
            aData.atomNumber=atoi(word.c_str());
             outfile<<word<<"  ";
             word=s.substr(13,2);
             aData.atomType=word.c_str();
             outfile<<word<<"  ";
             word=s.substr(17,3);
             aData.aminoAcidType=word.c_str();
            outfile<<word<<" ";
             word=s.substr(21,1);
             aData.chain=word.c_str();
             outfile<<word<<" ";
             word=s.substr(23,4);
             aData.residueNumber=atoi(word.c_str());
            outfile<<word<<"    ";
            word=s.substr(31,7);
             aData.xCoordinate=stod(word.c_str());
            outfile<<word<<" ";
            word=s.substr(39,7);
             aData.yCoordinate=stod(word.c_str());
            outfile<<word<<" ";
            word=s.substr(47,7);
             aData.zCoordinate=stod(word.c_str());
            outfile<<word<<"  ";
            word=s.substr(56,10);
            aData.unknown=word.c_str();
            outfile<<word<<"           ";
            word=s.substr(77,1);
            aData.unknown2=word.c_str();
            outfile<<word<<" ";
            atoms.push_back(aData);

        }

}

void extracthelixData(string s, ofstream & outfile)
{
    string word;
    helixData sData;

    int value;
     word=s.substr(0,5);
     string temp="";
        if(word=="HELIX")
        {
            temp=s.substr(19,1);

        }

         if(word=="HELIX"&& temp=="A")
        {
            outfile<<endl;
            sData.recordType=word;
            outfile<<word<<"   ";
            word=s.substr(8,2);
            outfile<<word<<" ";
            sData.serialnumber=atoi(word.c_str());
            word=s.substr(11,3);
           outfile<<word<<" ";
           sData.identifier=word;
           word=s.substr(15,3);
           outfile<<word<<" ";
           sData.residueName=word;
            word=s.substr(19,1);
            outfile<<word<<"  ";
            sData.chainIdentifier=word;
             word=s.substr(22,3);
             outfile<<word<<" ";
             sData.residueSequenceNumber=atoi(word.c_str());
             word=s.substr(26,4);
            outfile<<word<<" ";
            sData.unknown1=word;
             word=s.substr(31,1);
             outfile<<word<<" ";
             sData.terminalResidueName=word;
             word=s.substr(32,1);
             outfile<<word<<"";
             sData.chainidentifier=word;
             word=s.substr(34,3);
             outfile<<word<<" ";
             sData.residueSequenceNumber2=atoi(word.c_str());
             word=s.substr(38,38);
                outfile<<word<<" ";
             sData.unknown2=word;
             helix.push_back(sData);

        }

}

void extractsidechain(string s, ofstream & outfile)
{
    string word;
    atomData aData;
    int value;
     word=s.substr(0,5);
     string temp="";
         word=s.substr(0,4);
         string word1="";
         string word2="";
        if(word=="ATOM")
        {
            word1=s.substr(13,2);
        }

        if(((word=="ATOM") && (word1!="N "))&&((word=="ATOM") && (word1!="CA"))&&((word=="ATOM") && (word1!="C "))&&((word=="ATOM") && (word1!="O "))&&((word=="ATOM") && (word1!="H ")))
        {
            outfile<<endl;
            aData.atomName=word;
            outfile<<word<<"   ";
            word=s.substr(7,4);
            aData.atomNumber=atoi(word.c_str());
             outfile<<word<<"  ";
             word=s.substr(13,3);
             aData.atomType=word.c_str();
             outfile<<word<<" ";
             word=s.substr(17,3);
             aData.aminoAcidType=word.c_str();
            outfile<<word<<" ";
             word=s.substr(21,1);
             aData.chain=word.c_str();
             outfile<<word<<" ";
             word=s.substr(23,4);
             aData.residueNumber=atoi(word.c_str());
            outfile<<word<<"    ";
            word=s.substr(31,7);
             aData.xCoordinate=stod(word.c_str());
            outfile<<word<<" ";
            word=s.substr(39,7);
             aData.yCoordinate=stod(word.c_str());
            outfile<<word<<" ";
            word=s.substr(47,7);
             aData.zCoordinate=stod(word.c_str());
            outfile<<word<<"  ";
            word=s.substr(56,10);
            aData.unknown=word.c_str();
            outfile<<word<<"           ";
            word=s.substr(77,1);
            aData.unknown2=word.c_str();
            outfile<<word<<" ";
            atoms.push_back(aData);

        }

}

void extractsheetData(string s, ofstream & outfile)
{
  string word;
    sheetData sData;
    atomData aData;
    int value;
     word=s.substr(0,5);
     string temp="";
        if(word=="SHEET")
        {
            temp=s.substr(21,1);

        }

         if(word=="SHEET")
        {
           outfile<<endl;
            sData.recordType=word;
            outfile<<word<<" ";
            word=s.substr(17,3);
           outfile<<word<<" ";
           sData.residueName=word;
           word=s.substr(21,1);
           outfile<<word<<" ";
           sData.aminoType1=word;
            word=s.substr(22,4);
            outfile<<word<<" ";
            sData.sequenceNumber=atoi(word.c_str());
             word=s.substr(28,3);
             outfile<<word<<" ";
             sData.terminalResidueName=word;
             word=s.substr(32,1);
           outfile<<word<<" ";
           sData.aminoType2=word;
             word=s.substr(33,4);
             outfile<<word<<" ";
             sData.terminalISequenceNumber=atoi(word.c_str());

             sheets.push_back(sData);

        }
}
void extractBetastrand(string s, ofstream & outfile)
{
    string word;
    atomData aData;
    int value;
     word=s.substr(0,5);
     string temp="";
         word=s.substr(0,4);
         string word1="";
         string word2="";
        if(word=="ATOM")
        {
            word1=s.substr(13,2);
        }

        if(((word=="ATOM") && (word1=="N "))||((word=="ATOM") && (word1=="CA"))||((word=="ATOM") && (word1=="C ")))
        {
            outfile<<endl;
            aData.atomName=word;
            outfile<<word<<"   ";
            word=s.substr(7,4);
            aData.atomNumber=atoi(word.c_str());
             outfile<<word<<"  ";
             word=s.substr(13,3);
             aData.atomType=word.c_str();
             outfile<<word<<" ";
             word=s.substr(17,3);
             aData.aminoAcidType=word.c_str();
            outfile<<word<<" ";
             word=s.substr(21,1);
             aData.chain=word.c_str();
             outfile<<word<<" ";
             word=s.substr(23,4);
             aData.residueNumber=atoi(word.c_str());
            outfile<<word<<"    ";
            word=s.substr(31,7);
             aData.xCoordinate=stod(word.c_str());
            outfile<<word<<" ";
            word=s.substr(39,7);
             aData.yCoordinate=stod(word.c_str());
            outfile<<word<<" ";
            word=s.substr(47,7);
             aData.zCoordinate=stod(word.c_str());
            outfile<<word<<"  ";
            word=s.substr(56,10);
            aData.unknown=word.c_str();
            outfile<<word<<"           ";
            word=s.substr(77,1);
            aData.unknown2=word.c_str();
            outfile<<word<<" ";
            atoms.push_back(aData);

        }

}

