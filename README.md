### NOTE

ElConcatenero was discontinued in favor of TriFusion software

GitHub page: https://github.com/ODiogoSilva/TriFusion
Web page: http://odiogosilva.github.io/TriFusion/

#### About ElConcatenero3


ElConcatenero3 is a python program developed for the UNIX OS (but should work in any operating system with python 3.x) that converts and concatenates common population genetics and phylogenetics data file types (such as Fasta, Phylip and Nexus). The program works through command line exclusively, which makes it fast and scriptable, a usefull feature when dealing with lots of files simultaneously. It uses the ElParsito3.py module, which is where the file parsing and writting classes are.

It currently supports as input file formats:

- Fasta
- Phylip
- Nexus

And its able to convert/concatenate any of these file formats into:

- Fasta
- Phylip (Including the *.parts file required for multilocus datasets in RAxML)
- Nexus
- ZORRO weigths files (these weigths files that are associated with their respective alignment files can be jointly concatenated)

There is no need for instalation. The only requirement is that the ElParsito3.py module must be on the same directory as ElConcatenero3.py (or you can use any other way to let the main script know where the module is). I do recommend, to make it easier to call the program, that you add it to your $PATH variable. This can be done by declaring "export $PATH=$PATH:/path/to/ElConcatenero3" on your .bashrc, or your shell specific rc file.

Finally, please note that ElConcatero3.py is not immune to bugs, and I'll be happy to know about them through my e-mail (o.diogosilva@gmail.com) so that I can fix them. Any suggestions or comments are also welcome.

#### Options

ElConcatenero3.py has the following options (which can also be consulted by typing "ElConcatenero3.py -h" or "ElConcatenero3.py --help in the command line):

  -h, --help						**show this help message and exit**
 
  -g *GAP*							**Symbol for gap (default is '-')**
  
  -m *MISSING*						**Symbol for missing data (default is 'n')**
  
  -if *{fasta,nexus,phylip}*		**Format of the input file(s) (default is 'fasta')**
                        
  -of *{nexus,phylip,fasta}*		**Format of the ouput file (default is 'nexus')**
                        
  -c                    			**Used for convertion of the input files passed as**
									**arguments with the -in option. This flag precludes the**
									**usage of the -o option, as the output file name is**
									**automatically generated based on the input file name.**
                        
  -o *OUTFILE*            			**Name of the output file**
  
  -in *INFILE [INFILE ...]*			**Provide the input file name. If multiple files are**
									**provided, plase separated the names with spaces**
  -z								**Use this option if you wish to concatenate auxiliary**
									**Zorro files associated with each alignment. Note that**
									**the auxiliary files must have the same prefix of the**
									**alignment file, with the addition of '_zorro.out'**
								
  -zfile *[ZORRO_INFILE ...]*		** provide the sufix for the concatenated zorro file**
									**(default is '_zorro.out')**

  -code *{DNA, Protein}*			**The coding of your molecular sequences**
  
  -rm *[REMOVE ...]*				**Removes the specified taxa from the final alignment.**
									**Multiple taxa mais be specified and separated by**
									**whitespace**
								
  -\-pickle-taxa *{dump,load}*		**dump option: Only output a picke object with the taxa**
									**names of the input alignment; load option: loads the**
									**taxa names from a pickle object to be incorporated in**
									**the output alignment**
									
  -\-check 							**Provides a final check for the lengths of the**
									**alignment sequences**

#####Note: The order of the options does not matter.
		
#### Usage

##### Conversion (Fasta to Nexus)

ElConcatenero3.py -c -of nexus -in input_file.fas

##### Conversion (Phylip to Fasta)

ElConcatenero3.py -c -if phylip -of nexus -in input_file.phy

##### Conversion (Nexus to Phylip and fasta)

ElConcatenero3.py -c -if nexus -of phylip fasta -in input_file.nex

##### Conversion (Multiple Fasta files to Nexus)

ElConcatenero3.py -c -of nexus -in input_file1.fas input_file2.fas input_file3.fas input_file4.fas [...] input_fileN.fas

or

ElConcatenero3.py -c -of nexus -in *.fas

##### Concatenation (Nexus to Nexus)

ElConcetenero.py -if nexus -of nexus -in input_file1.nex input_file2.nex input_file3.nex (...) input_fileN.nex -o concatenated_file

or

ElConcetenero.py -if nexus -of nexus -in *.nex -o concatenated_file

##### Concatenation (Fasta to Phylip)

ElConcatenero3.py -if fasta -of phylip -in input_file1.fas input_file2.fas (...) input_fileN.fas -o concatenated_file

or

ElConcatenero3.py -if fasta -of phylip -in *.fas -o concatenated_file

#### ToDo

- Ability to parse concatenated datasets, either to separate them or to further concatenate
- Auto recognition of DNA or Protein sequences
- Auto recognition of input format
