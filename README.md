#### About ElConcatenero

ElConcatenero is a python program developed for UNIX OS that converts and concatenates common population genetics and phylogenetics data file types (such as FASTA, Phylip, Nexus, IMa2). The program works through command line exclusively, which makes it fast and scriptable, a usefull feature when dealing with lots of files simultaneously as I do. It uses the ElParsito.py module, which is where the file parsing functions are.

It currently supports as input file formats:

- FASTA
- Phylip
- Nexus

And its able to convert/concatenate any of these file formats into:

- Fasta
- Phylip (Including the *.parts file required for multilocus datasets in RAxML)
- Nexus
- IMa2 input format

There is no need for instalation. The only requirement is that the ElParsito.py module must be on the same directory as ElConcatenero.py (or you can use any other way to let the main script know where the module is). I do recommend, to make it easier to call the program, that you add it to your $PATH variable. 

Finally, please note that ElConcatero.py is far from immune to bugs and crashes, and I'll be happy to know about them through my e-mail (o.diogosilva@gmail.com) so that I can fix them. Any suggestions or comments are also welcome.

#### Options

ElConcatenero.py has the following options (which can also be consulted by typing "ElConcatenero.py -h" in the command line):

  -h, --help					**show this help message and exit**
 
  -g *GAP*						**Symbol for gap (default is '-')**
  
  -m *MISSING*					**Symbol for missing data (default is 'n')**
  
  -if *{fasta,nexus,phylip}*		**Format of the input file(s) (default is 'fasta')**
                        
  -of *{nexus,phylip,fasta,ima2}*	**Format of the ouput file (default is 'nexus')**
                        
  -c                    		**Used for convertion of the input files passed as**
								**arguments with the -in option. This flag precludes the**
								**usage of the -o option, as the output file name is**
								**automatically generated based on the input file name.**
                        
  -o *OUTFILE*            		**Name of the output file**
  
  -in *INFILE [INFILE ...]*		**Provide the input file name. If multiple files are**
								**provided, plase separated the names with spaces**
								
#####Note:

The order of the options does not matter.
								
#### Usage

##### Conversion (FASTA to Nexus)

ElConcatenero.py -c -if fasta -of nexus -in input_file.fas

##### Conversion (Phylip to FASTA)

ElConcatenero.py -c if phylip -of nexus -in input_file.phy

##### Concatenation (Nexus to Nexus)

ElConcetenero.py -if nexus -of nexus -in input_file1.nex input_file2.nex input_file3.nex (...) input_fileN.nex -o concatenated_file

##### Concatenation (FASTA to Phylip)

ElConcatenero.py -if fasta -of phylip -in input_file1.fas input_file2.fas (...) input_fileN.fas -o concatenated_file

##### Concatenation (Phylip to IMa2) - The user is prompted for additional information, required to produce the output file in IMa2 format

ElConcatenero.py -if phylip -of ima2 -in input_file1.phy input_file2.phy (...) input_fileN.phy -o Ima2_file

#### Tips on using ElConcatenero with multiple files

These tips only apply to UNIX OS, as they resort to bash.

- Assuming you are in a directory with N FASTA files, which you would like to convert to Nexus, you could do something like:

[starting pseudo-code]

for file in `ls | grep .fas*`;

 do path/to/ElConcatenero.py -c -if fasta -of nexus -in $file;
 
done;

[ending pseudo-code]

And you'll have a potentially huge ammount of files converted in no time.

- You may also have several directories each containing multiple nexus files that you wish to concatenate separately. You could do somenting like (assuming you are on the parent directory):

[starting pseudo-code]
for directory in `ls -d */`;

 do cd $directory
 
 /path/to/ElConcatenero.py -if nexus -of nexus -in *.nex -o concatenated_dataset
 
 cd ../
 
done

[ending pseudo-code]

#### ToDo

- Ability to parse concatenated datasets, either to separate them or to further concatenate
