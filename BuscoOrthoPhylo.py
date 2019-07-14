#!/usr/bin/python3

'''
Written by Stephen A. Wyka (2019)

This is a wrapper script that takes multiple BUSCO output folders or the output folder 
from Orthofinder to concatenate, align, and contruct phylogentic tree(s) of single copy 
BUSCOs/orthologs. Please be aware of directory strucuture required for BUSCO innput.

For BUSCO use please set up your -d/--directory input like so:

Parent_directory
    |
    |
    |
    +--Genus_species  <- This name will be used in FASTA headers and Phylogenetic trees
    |       |
    |       |
    |       +--run_xxx  <- This is your BUSCO output folder
    |
    |
    +--Genus_species
            |
            |
            +--run_xxx

For Orthofinder use:

Parent_directory = Results_XXX output folder from Orthofinder

'''

import os, sys, re, argparse, shutil, textwrap, subprocess, random
from collections import defaultdict, OrderedDict
currentdir = os.getcwd()
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

class MyFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(
    usage='./%(prog)s [options] -d directory -o output_directory -i busco/orthofinder',
    description = '''    Wrapper to take multiple BUSCO output folders or output folder 
    from Orthofinder and concatenates, aligns, and constructs phylogenetic tree(s) of 
    single copy BUSCOs/orthologs.''',
    
    epilog = """Written by Stephen A. Wyka (2019)""",
    formatter_class = MyFormatter)

parser.add_argument(
    '-d',
    '--directory',
    required=True,
    help = 'Parent directory containing directory for each species',
    metavar=''
)
parser.add_argument(
    '-i',
    '--input',
    required=True,
    choices=['busco','orthofinder'],
    help = 'Program used for input directory data [busco|orthofinder]',
    metavar=''
)
parser.add_argument(
    '-o',
    '--out',
    required=True,
    help = 'Name for output directory',
    metavar=''
)
parser.add_argument(
    '-k',
    '--key',
    help = 'Keyword to help find correct BUSCO output folder [ex. -k sordario when '\
    'busco output folder is run_sordario',
    metavar=''
)
parser.add_argument(
    '-t',
    '--type',
    default='prot',
    choices = ['prot', 'nuc', 'both'],
    help = 'Choice of sequence output [prot|nuc|both] (only prot available when '\
    'using orthofinder)[default: prot]',
    metavar=''
)
parser.add_argument(
    '-sub',
    '--subsample',
    type=float,
    help = 'Percent of genes to subsample for a preliminary analysis [ex. 0.10]',
    metavar=''
)
parser.add_argument(
    '-l',
    '--list',
    help = 'File contains list of single copy BUSCOs or Orthogroups for Phylogeny'\
    ' on (one per line)',
    metavar=''
)
parser.add_argument(
    '--stop',
    choices = ['ac', 'aa'],
    help = 'Stop after concatenation [ac] stop after alignment [aa]',
    metavar=''
)
parser.add_argument(
    '--resume',
    choices = ['ac', 'aa'],
    help = 'Resume after concatenation [ac] resume after alignment [aa]',
    metavar=''
)
parser.add_argument(
    '-c',
    '--cpus',
    default=2,
    type=int,
    help = 'Cores to use [default: 2]',
    metavar=''
)
parser.add_argument(
    '-a',
    '--aligner',
    default='mafft',
    choices = ['mafft', 'muscle'],
    help = 'Choice of sequence aligner [mafft|muscle] [default: mafft]',
    metavar=''
)
parser.add_argument(
    '-ma_a',
    '--mafft_args',
    nargs='+',
    help = "Additional flags, separated by commas, within quotes, no spaces after commas, "\
    "to be used in MAFFT alignment [ex. 'localpair,maxiterate 0-1000']",
    metavar=''
)
parser.add_argument(
    '-mu_i',
    '--muscle_iter',
    type=int,
    default=16,
    help = 'Maximum number of iterations for MUSCLE [default: 16]',
    metavar=''
)
parser.add_argument(
    '-mu_d',
    '--muscle_diags',
    action='store_true',
    help = 'Turn on find diagonals for MUSCLE (faster for similar sequences) [default: OFF]',
)
parser.add_argument(
    '-ph',
    '--phylogeny',
    default='fasttree',
    choices = ['fasttree', 'raxml'],
    help = 'Choice of sequence output [fasttree|raxml] [default: fasttree]',
    metavar=''
)
parser.add_argument(
    '-b',
    '--boot',
    type=int,
    default=1000,
    help = 'Number of bootstraps. Seriously consider reducing if using raxml. [default: 1000]',
    metavar=''
)
parser.add_argument(
    '-ft_a',
    '--fasttree_args',
    nargs='+',
    help = "Additional flags, separated by commas, within quotes, no spaces after commas,"\
    "to be used in fasttree run [ex. 'slow,bionj,spr 4']",
    metavar=''
)
parser.add_argument(
    '-rx_a',
    '--raxml_args',
    nargs='+',
    help = "Additional flags, separated by commas, within quotes, no spaces after commas,"\
    "to be used in raxml run [ex. 'k,f a,t file.nwk']",
    metavar=''
)
parser.add_argument(
    '-rx_og',
    '--raxml_outgroup',
    nargs='+',
    help = "Names of outgroup species to use in raxml run [ex. 'rat,mouse,fly,human']. "\
    "Names must correlate to headers in alignment file",
    metavar=''
)
parser.add_argument(
    '-rx_p_sub',
    '--raxml_prot_sub',
    default='PROTCATAUTO',
    help = 'Protein substitution model to use [see raxmlHPC -h for all -m flag options] '\
    '[default: PROTCATAUTO]',
    metavar=''
)
parser.add_argument(
    '-rx_n_sub',
    '--raxml_nuc_sub',
    default='GTRCAT',
    help = 'Nucleotide substitution model to use [see raxmlHPC -h for all -m flag options] '\
    '[default: GTRCAT]',
    metavar=''
)
parser.add_argument(
    '-rx_al',
    '--raxml_algorithm',
    default='a',
    help = 'Algorithm to use [see raxmlHPC -h for all -f flag options] '\
    '[default: a]',
    metavar=''
)
parser.add_argument(
    '-rx_t',
    '--raxml_threads',
    action='store_true',
    help = 'Turn off PTHREADS [default: ON]',
)
args = parser.parse_args()

# Load functions

def which_path(file_name):
    for path in os.environ["PATH"].split(os.pathsep):
        full_path = os.path.join(path, file_name)
        if os.path.exists(full_path) and os.access(full_path, os.X_OK):
            return full_path
    return None

def walklevel(input_dir, level=1):
    input_dir = input_dir.rstrip(os.path.sep)
    if os.path.isdir(input_dir):
        assert os.path.isdir(input_dir)
        num_sep = input_dir.count(os.path.sep)
        for root, dirs, files in os.walk(input_dir):
            yield root, dirs, files
            num_sep_this = root.count(os.path.sep)
            if num_sep + level <= num_sep_this:
                del dirs[:]
    else:
        None

def busco_output(search_dir):
    try:
        if not args.key:
            busco_folder_count = 0
            for dir in os.listdir(search_dir):
                if dir.startswith('run_'):
                    busco_folder_count += 1
                    busco_out = os.path.abspath(os.path.join(search_dir, dir))
                    if busco_folder_count > 1:
                        raise 
    except:
        print('ERROR: Multiple BUSCO output folders found per species folder, please',\
        'into a keyword with -k or --key to help select which BUSCO output folder to use')
        sys.exit(1)
    else:
        try:
            if args.key:
                for dir in os.listdir(search_dir):
                    if dir.startswith('run_') and args.key in dir:
                        busco_out = os.path.abspath(os.path.join(search_dir, dir))
            return busco_out
        except:
            print('ERROR: Problem with BUSCO output folder names. Keyword is not is some of '\
            'the folders')
            sys.exit(1)

# Create output folder and define paths

if not os.path.isdir(args.out):
    os.makedirs(os.path.join(args.out, 'species'))
    os.makedirs(os.path.join(args.out, 'alignment'))
    os.makedirs(os.path.join(args.out, 'phylogeny'))
result_dir = os.path.abspath(os.path.join(currentdir, args.out))
species_dir = os.path.abspath(os.path.join(result_dir, 'species'))
align_dir = os.path.abspath(os.path.join(result_dir, 'alignment'))
phylo_dir = os.path.abspath(os.path.join(result_dir, 'phylogeny'))
if not os.path.isdir(species_dir):
    dirs = [os.path.join(result_dir, 'species'), os.path.join(result_dir, 'alignment'),
    os.path.join(result_dir, 'phylogeny'), os.path.join(currentdir, args.out)]
    for d in dirs:
        if not os.path.isdir(d):
            os.makedirs(d)
input_dir = os.path.abspath(os.path.join(currentdir, args.directory))

# Checking dependencies

if args.aligner == 'muscle':
    try:
        if which_path('muscle'):
            MUSCLE = 'muscle'
        else:
            raise
    except:
        print('ERROR: muscle not found, please make sure parent directory of',\
        'muscle is located in $PATH')
        sys.exit(1)

if args.aligner == 'mafft':
    try:
        if which_path('mafft'):
            MAFFT = 'mafft'
        else:
            raise
    except:
        print('ERROR: mafft not found, please make sure parent directory of',\
        'mafft is located in $PATH')
        sys.exit(1)

if args.phylogeny == 'fasttree':
    try:
        if which_path('FastTree'):
            FASTTREE = 'FastTree'
        else:
            raise
    except:
        print('ERROR: FastTree not found, please make sure parent directory of',\
        'FastTree is located in $PATH')
        sys.exit(1)

if args.phylogeny == 'raxml':
    try:
        if which_path('raxmlHPC'):
            RAXML = 'raxmlHPC'
        else:
            raise
    except:
        print('ERROR: raxmlHPC not found, please make sure parent directory of',\
        'raxmlHPC is located in $PATH')
        sys.exit(1)

# Creating paths for some arguments

if args.list:
    filter_list = os.path.abspath(os.path.join(currentdir,args.list))

# Scripts for BUSCO input
if not args.resume: # If resume is called, this will skip concatenation step.
    if args.input == 'busco':

        input_directory = [os.path.abspath(x[0]) for x in walklevel(input_dir)]
        input_directory.remove(os.path.abspath(input_dir))
        for dir in os.listdir(species_dir):
            shutil.rmtree(os.path.join(species_dir,dir))
        for i in input_directory:
            species_name = os.path.basename(i)
            if not os.path.isdir(os.path.join(species_dir, species_name)):
                os.makedirs(os.path.join(species_dir, species_name))

    # Create list of BUSCOS

        print('Creating dictionary of BUSCOs per species')

        buscos = OrderedDict()
        for i in input_directory:
            species_name = os.path.basename(i)
            buscos[species_name] = []
            busco_out = busco_output(i)
            for path, directories, files in os.walk(os.path.join(busco_out, 'single_copy_busco_sequences')):
                for file in files:
                    fpath = os.path.join(path, file)
                    name, extension = file.split('.')
                    buscos[species_name] = buscos.get(species_name, []) + [name]

    # Finding common BUSCOs

        print('Finding common BUSCOs between all species')
        common_buscos = []
        lengths = []
        for values in buscos.values():
            lengths.append(len(values))

        minLength = min(lengths)

        for keys, values in buscos.items():
            if int(len(values)) == int(minLength):
                common_buscos = set(sorted(values))
            else:
                continue

        for keys, values in buscos.items():
            common_buscos = set(sorted(common_buscos)) & set(sorted(values))
        

    # Turn common BUSOCs into dictionary for faster searching and get subsample if applied

        common_dict = OrderedDict()
        if args.subsample:
            sub_number = round(len(common_buscos) * args.subsample)
            sub_busco = random.sample(common_buscos, sub_number)
            print('Subsampling {} of {} shared BUSCOs for downstream analysis'.format(sub_number, len(common_buscos)))
            for x in sub_busco:
                common_dict[x] = []
        elif args.list:
            with open(filter_list, 'r') as filter:
                given_list = [x.strip() for x in filter]
                try:
                    test_list = list(set(given_list)-set(common_buscos))
                    if test_list:
                        raise
                except:
                    print("ERROR: Provided list contains BUSCOs for which provided species do not"\
                    " share in common. Problematic BUSCOs are:")
                    print(','.join(test_list))
                    sys.exit(1)
                for x in given_list:
                    common_dict[x] = []
                print('Using {} provided BUSCOs shared in all {}species for downstream'\
                'analysis'.format(len(given_list),len(buscos.keys())))
        else:
            for x in common_buscos:
                common_dict[x] = []
            print('Using {} BUSCOs shared in all {} species for downstream'\
            'analysis'.format(len(common_buscos),len(buscos.keys())))

    # Creating FASTA files for each species

        if args.type == 'prot':
            print('Creating protein FASTA files of BUSCOs for each species')
            for i in input_directory:
                species_name = os.path.basename(i)
                species_results = os.path.abspath(os.path.join(species_dir,species_name))
                single_prot = os.path.abspath(os.path.join(species_results,species_name+'_single_busco_prot.fasta'))
                concat_prot = os.path.abspath(os.path.join(species_results,species_name+'_concat_busco_prot.fasta'))
                busco_out = busco_output(i)
                with open(single_prot, 'w') as s_prot_out, open(concat_prot, 'w') as c_prot_out:
                    c_prot_out.write('>'+species_name + '\n')
                    concat_list = OrderedDict()
                    for path, directories, files in os.walk(os.path.join(busco_out, 'single_copy_busco_sequences')):
                        for file in files:
                            fpath = os.path.join(path, file)
                            name, extension = file.split('.')
                            if name in common_dict.keys() and extension == 'faa':
                                with open(fpath, 'r') as prot_file:
                                    prots = prot_file.readlines()
                                    s_prot_out.write(''.join(prots))
                                    concat_list[prots[1].strip()] = []
                    c_prot_out.write(textwrap.fill(''.join(list(concat_list.keys())),width=80))
        if args.type == 'nuc':
            print('Creating nucleotide FASTA files of BUSCOs for each species')
            for i in input_directory:
                species_name = os.path.basename(i)
                species_results = os.path.abspath(os.path.join(species_dir,species_name))
                single_nuc = os.path.abspath(os.path.join(species_results,species_name+'_single_busco_nuc.fasta'))
                concat_nuc = os.path.abspath(os.path.join(species_results,species_name+'_concat_busco_nuc.fasta'))
                busco_out = busco_output(i)
                with open(single_nuc, 'w') as s_nuc_out, open(concat_nuc, 'w') as c_nuc_out:
                    c_nuc_out.write('>'+species_name + '\n')
                    concat_list = OrderedDict()
                    for path, directories, files in os.walk(os.path.join(busco_out, 'single_copy_busco_sequences')):
                        for file in files:
                            fpath = os.path.join(path, file)
                            name, extension = file.split('.')
                            if name in common_dict.keys() and extension == 'fna':
                                with open(fpath, 'r') as nuc_file:
                                    nucs = nuc_file.readlines()
                                    s_nuc_out.write(nucs[0] + nucs[1].upper())
                                    concat_list[nucs[1].strip().upper()] = []
                    c_nuc_out.write(textwrap.fill(''.join(list(concat_list.keys())),width=80))
        if args.type =='both':
            print('Creating protein and nucleotide FASTA files of BUSCOs for each species')
            for i in input_directory:
                species_name = os.path.basename(i)
                species_results = os.path.abspath(os.path.join(species_dir,species_name))
                single_prot = os.path.abspath(os.path.join(species_results,species_name+'_single_busco_prot.fasta'))
                concat_prot = os.path.abspath(os.path.join(species_results,species_name+'_concat_busco_prot.fasta'))
                single_nuc = os.path.abspath(os.path.join(species_results,species_name+'_single_busco_nuc.fasta'))
                concat_nuc = os.path.abspath(os.path.join(species_results,species_name+'_concat_busco_nuc.fasta'))
                busco_out = busco_output(i)
                with open(single_prot, 'w') as s_prot_out, open(concat_prot, 'w') as c_prot_out:
                    c_prot_out.write('>'+species_name + '\n')
                    concat_list = OrderedDict()
                    for path, directories, files in os.walk(os.path.join(busco_out, 'single_copy_busco_sequences')):
                        for file in files:
                            fpath = os.path.join(path, file)
                            name, extension = file.split('.')
                            if name in common_dict.keys() and extension == 'faa':
                                with open(fpath, 'r') as prot_file:
                                    prots = prot_file.readlines()
                                    s_prot_out.write(''.join(prots))
                                    concat_list[prots[1].strip()] = []
                    c_prot_out.write(textwrap.fill(''.join(list(concat_list.keys())),width=80))
                with open(single_nuc, 'w') as s_nuc_out, open(concat_nuc, 'w') as c_nuc_out:
                    c_nuc_out.write('>'+species_name + '\n')
                    concat_list = OrderedDict()
                    for path, directories, files in os.walk(os.path.join(busco_out, 'single_copy_busco_sequences')):
                        for file in files:
                            fpath = os.path.join(path, file)
                            name, extension = file.split('.')
                            if name in common_dict.keys() and extension == 'fna':
                                with open(fpath, 'r') as nuc_file:
                                    nucs = nuc_file.readlines()
                                    s_nuc_out.write(nucs[0] + nucs[1].upper())
                                    concat_list[nucs[1].strip().upper()] = []
                    c_nuc_out.write(textwrap.fill(''.join(list(concat_list.keys())),width=80))

    # Creating final concatenation file

        final_prot = os.path.abspath(os.path.join(align_dir,'Concatenated_{}_proteins.fasta'.format(args.input)))
        final_nuc = os.path.abspath(os.path.join(align_dir,'Concatenated_{}_nucleotides.fasta'.format(args.input)))
        if args.type == 'prot':
            print('Creating final concatenated protein FASTA file')
            with open(final_prot, 'w') as f_prot:
                for i in input_directory:
                    species_name=os.path.basename(i)
                    species_results = os.path.abspath(os.path.join(species_dir,species_name))
                    single_prot = os.path.abspath(os.path.join(species_results,species_name+'_single_busco_prot.fasta'))
                    concat_prot = os.path.abspath(os.path.join(species_results,species_name+'_concat_busco_prot.fasta'))
                    with open(concat_prot, 'r') as c_prot:
                        spec_prot = c_prot.read()
                        f_prot.write(spec_prot + '\n')
        if args.type == 'nuc':
            print('Creating final concatenated nucleotide FASTA file')
            with open(final_nuc, 'w') as f_nuc:
                for i in input_directory:
                    species_name=os.path.basename(i)
                    species_results = os.path.abspath(os.path.join(species_dir,species_name))
                    single_nuc = os.path.abspath(os.path.join(species_results,species_name+'_single_busco_nuc.fasta'))
                    concat_nuc = os.path.abspath(os.path.join(species_results,species_name+'_concat_busco_nuc.fasta'))
                    with open(concat_nuc, 'r') as c_nuc:
                        spec_nuc = c_nuc.read()
                        f_nuc.write(spec_nuc + '\n')
        if args.type =='both':
            print('Creating final concatenated protein and nucleotide FASTA files')
            with open(final_prot, 'w') as f_prot, open(final_nuc, 'w') as f_nuc:
                for i in input_directory:
                    species_name=os.path.basename(i)
                    species_results = os.path.abspath(os.path.join(species_dir,species_name))
                    single_prot = os.path.abspath(os.path.join(species_results,species_name+'_single_busco_prot.fasta'))
                    concat_prot = os.path.abspath(os.path.join(species_results,species_name+'_concat_busco_prot.fasta'))
                    single_nuc = os.path.abspath(os.path.join(species_results,species_name+'_single_busco_nuc.fasta'))
                    concat_nuc = os.path.abspath(os.path.join(species_results,species_name+'_concat_busco_nuc.fasta'))
                    with open(concat_prot, 'r') as c_prot:
                        spec_prot = c_prot.read()
                        f_prot.write(spec_prot + '\n')
                    with open(concat_nuc, 'r') as c_nuc:
                        spec_nuc = c_nuc.read()
                        f_nuc.write(spec_nuc + '\n')

    # Scripts for Orthofinder input 
    ## Would be nice to reduce code and combine some ascept (for future thought)

    if args.input == 'orthofinder':

        if args.type == 'nuc' or args.type == 'both':
            print('ERROR: We cannot run nucleotide alignment and phylogeny on orthofinder data',\
            'as orthofinder does not provide nucleotide FASTA files')
            sys.exit(1)
        common_ortho = []
        common_ortho_dict = OrderedDict()
        species_dict = defaultdict(str)
        species_ortho_dict = defaultdict(str)
        version = []
        if os.path.isdir(os.path.join(input_dir,'Single_Copy_Orthologue_Sequences')): 
            version = 2.3
        else:
            version = 2.2

    # Get list of common ortholog names
        
        print('Creating dictionary of single copy orthologs')
        if version == 2.3:
            single_seqs = os.path.abspath(os.path.join(input_dir, 'Single_Copy_Orthologue_Sequences'))
            for path, dirs, files in os.walk(single_seqs):
                for file in files:
                    name, extension = file.split('.')
                    common_ortho.append(name)
        if version == 2.2:
            single_copy = os.path.abspath(os.path.join(input_dir, 'SingleCopyOrthogroups.txt'))
            with open(single_copy, 'r') as s_ortho:
                for ortho in s_ortho:
                    common_ortho.append(ortho.strip())

    # Create ortholog dictionary and subsample if called

        if args.subsample:
            sub_number = round(len(common_ortho) * args.subsample)
            sub_ortho = random.sample(list(common_ortho), sub_number)
            print('Subsampling {} of {} shared Orthologs for downstream analysis'.format(sub_number, len(common_ortho)))
            for x in sub_ortho:
                common_ortho_dict[x] = []
        elif args.list:
            with open(filter_list, 'r') as filter:
                given_list = [x.strip() for x in filter]
                try:
                    test_list = list(set(given_list)-set(common_ortho))
                    if test_list:
                        raise
                except:
                    print("ERROR: Provided list contains ortho groups for which provided species do not"\
                    " share in common. Problematic ortho groups are:") 
                    print(','.join(test_list))
                    sys.exit(1)
                for x in given_list:
                    common_ortho_dict[x] = []
                print('Using {} provided orthologs shared in all species for downstream analysis'.format(len(given_list)))
        else:
            for x in common_ortho:
                common_ortho_dict[x] = []
            print('Using {} Orthologs shared in all species for downstream analysis'.format(len(common_ortho)))

        # if version == 2.3:
             # do nothing for now
        if version == 2.2:
            for path, dirs, files in os.walk(input_dir):
                if os.path.isdir(os.path.join(path,'Sequences')):
                    single_seqs = os.path.abspath(os.path.join(path,'Sequences'))
                    for path, dirs, files in os.walk(single_seqs):
                        for file in files:
                            fpath = os.path.join(path, file)
                            name, extension = file.split('.')
                            if name in common_ortho_dict.keys():
                                with open(fpath, 'r') as seq_file:
                                    name = ''
                                    for line in seq_file:
                                        if line.startswith('>'):
                                            ortholog=line.strip().replace('>','')
                                            name=line.split('_')[0].replace('>','')
                                            continue
                                        species_dict[name]+=line.strip()
                                        species_ortho_dict[ortholog]+=line.strip()

    # Create folders for species found 

        for dir in os.listdir(species_dir):
            shutil.rmtree(os.path.join(species_dir,dir))
        for i in species_dict.keys():
            if not os.path.isdir(os.path.join(species_dir, i)):
                os.makedirs(os.path.join(species_dir, i))

    # Write FASTA files for individual species and concatenations
        
        print('Creating FASTA files of orthologs for each species and final concatenations')
        final_prot = os.path.abspath(os.path.join(align_dir,'Concatenated_{}_proteins.fasta'.format(args.input)))
        with open(final_prot, 'w') as cf_prot:
            for key, value in species_dict.items():
                cf_prot.write('>' + key + '\n' + textwrap.fill(value,width=80) + '\n')
        input_directory = [os.path.abspath(x[0]) for x in walklevel(species_dir)]
        input_directory.remove(os.path.abspath(species_dir))
        for i in input_directory:
            species_name=os.path.basename(i)
            species_results = os.path.abspath(os.path.join(species_dir, species_name))
            single_prot = os.path.abspath(os.path.join(species_results, species_name+'_single_orthos.fasta'))
            concat_prot = os.path.abspath(os.path.join(species_results, species_name+'_concat_orthos.fasta'))
            with open(concat_prot, 'w') as c_prot, open(single_prot, 'w') as s_prot:
                for key, value in species_dict.items():
                    if key == species_name:
                        c_prot.write('>' + key + '\n' + textwrap.fill(value,width=80) + '\n')
                for key, value in species_ortho_dict.items():
                    if species_name == key.split('_')[0]:
                        s_prot.write('>' + key + '\n' + textwrap.fill(value,width=80) + '\n')

# Stop after concatenation if called

if args.stop == 'ac':
    print('Done')
    sys.exit(1)

# Run alignment program
if not args.resume == 'aa': # If resume is called, but not aa we will run aligner
    print('Running {} aligner on concatenated FASTA file(s)'.format(args.aligner))
    for file in os.listdir(align_dir):
        if args.input == 'busco':
            if not file.startswith('Concatenated_busco'):
                os.remove(os.path.join(align_dir,file))
        if args.input == 'orthofinder':
            if not file.startswith('Concatenated_orthofinder'):
                os.remove(os.path.join(align_dir,file))
    final_prot = os.path.abspath(os.path.join(align_dir,'Concatenated_{}_proteins.fasta'.format(args.input)))
    final_nuc = os.path.abspath(os.path.join(align_dir,'Concatenated_{}_nucleotides.fasta'.format(args.input)))
    prot_align = os.path.abspath(os.path.join(align_dir,'{}_{}_prot_alignment.fasta'.format(args.aligner, args.input)))
    nuc_align = os.path.abspath(os.path.join(align_dir,'{}_{}_nuc_alignment.fasta'.format(args.aligner, args.input)))
    muscle_prot_log = os.path.abspath(os.path.join(align_dir,'muscle_{}_prot.log'.format(args.input)))
    muscle_nuc_log = os.path.abspath(os.path.join(align_dir,'muscle_{}_prot.log'.format(args.input)))
    if args.type == 'prot':
        if args.aligner == 'mafft':
            if args.mafft_args:
                mafft_args = ''.join(args.mafft_args).split(',')
                mafft_args = [a.replace(a,'--'+a) for a in mafft_args]
                cmd = [MAFFT,'--thread',str(args.cpus),' '.join(mafft_args),final_prot,'>',prot_align]
                os.system(' '.join(cmd))
            else:
                cmd = [MAFFT,'--thread',str(args.cpus),final_prot,'>',prot_align]
                os.system(' '.join(cmd))
        if args.aligner == 'muscle':
            if args.muscle_diags:
                cmd = [MUSCLE,'-diags','-maxiters',str(args.muscle_iter),'-log',muscle_prot_log,
                '-in',final_prot,'-out',prot_align]
                os.system(' '.join(cmd))
            else:
                cmd = [MUSCLE,'-maxiters',str(args.muscle_iter),'-log',muscle_prot_log,
                '-in',final_prot,'-out',prot_align]
                os.system(' '.join(cmd))
    if args.type == 'nuc':
        if args.aligner == 'mafft':
            if args.mafft_args:
                mafft_args = ''.join(args.mafft_args).split(',')
                mafft_args = [a.replace(a,'--'+a) for a in mafft_args]
                cmd = [MAFFT,'--thread',str(args.cpus),' '.join(mafft_args),final_nuc,'>',nuc_align]
                os.system(' '.join(cmd))
            else:
                cmd = [MAFFT,'--thread',str(args.cpus),final_nuc,'>',nuc_align]
                os.system(' '.join(cmd))
        if args.aligner == 'muscle':
            if args.muscle_diags:
                cmd = [MUSCLE,'-diags','-maxiters',str(args.muscle_iter),'-log',muscle_nuc_log,
                '-in',final_nuc,'-out',nuc_align]
                os.system(' '.join(cmd))
            else:
                cmd = [MUSCLE,'-maxiters',str(args.muscle_iter),'-log',muscle_nuc_log,
                '-in',final_nuc,'-out',nuc_align]
                os.system(' '.join(cmd))
    if args.type == 'both':
        if args.aligner == 'mafft':
            if args.mafft_args:
                mafft_args = ''.join(args.mafft_args).split(',')
                mafft_args = [a.replace(a,'--'+a) for a in mafft_args]
                cmd = [MAFFT,'--thread',str(args.cpus),' '.join(mafft_args),final_prot,'>',prot_align]
                os.system(' '.join(cmd))
                cmd = [MAFFT,'--thread',str(args.cpus),' '.join(mafft_args),final_nuc,'>',nuc_align]
                os.system(' '.join(cmd))
            else:
                cmd = [MAFFT,'--thread',str(args.cpus),final_prot,'>',prot_align]
                os.system(' '.join(cmd))
                cmd = [MAFFT,'--thread',str(args.cpus),final_nuc,'>',nuc_align]
                os.system(' '.join(cmd))
        if args.aligner == 'muscle':
            if args.muscle_diags:
                cmd = [MUSCLE,'-diags','-maxiters',str(args.muscle_iter),'-log',muscle_prot_log,
                '-in',final_prot,'-out',prot_align]
                os.system(' '.join(cmd))
                cmd = [MUSCLE,'-diags','-maxiters',str(args.muscle_iter),'-log',muscle_nuc_log,
                '-in',final_nuc,'-out',nuc_align]
                os.system(' '.join(cmd))
            else:
                cmd = [MUSCLE,'-maxiters',str(args.muscle_iter),'-log',muscle_prot_log,
                '-in',final_prot,'-out',prot_align]
                os.system(' '.join(cmd))
                cmd = [MUSCLE,'-maxiters',str(args.muscle_iter),'-log',muscle_nuc_log,
                '-in',final_nuc,'-out',nuc_align]
                os.system(' '.join(cmd))

# Stop after alignments if called

if args.stop == 'aa':
    print('Done')
    sys.exit(1)

# Run phylogeny

print('Running {} on {} alignment file(s)'.format(args.phylogeny,args.aligner))
for file in os.listdir(phylo_dir):
    os.remove(os.path.join(phylo_dir,file))
prot_tree = os.path.abspath(os.path.join(phylo_dir,'{}_{}_{}_prot_tree.nwk'.format(
    args.phylogeny,args.aligner,args.input)))
nuc_tree = os.path.abspath(os.path.join(phylo_dir,'{}_{}_{}_nuc_tree.nwk'.format(
    args.phylogeny,args.aligner,args.input)))
prot_align = os.path.abspath(os.path.join(align_dir,'{}_{}_prot_alignment.fasta'.format(args.aligner, args.input)))
nuc_align = os.path.abspath(os.path.join(align_dir,'{}_{}_nuc_alignment.fasta'.format(args.aligner, args.input)))
if args.phylogeny == 'fasttree':
    fasttree_log = os.path.abspath(os.path.join(phylo_dir,'{}_{}_fasttree.log'.format(args.aligner, args.input)))
    if args.type == 'prot':
        if args.fasttree_args:
            fast_args = ''.join(args.fasttree_args).split(',')
            fast_args = [a.replace(a,'-'+a) for a in fast_args]
            cmd = [FASTTREE,'-seed','12345','-gamma','-boot',str(args.boot),' '.join(fast_args),'-log',
            fasttree_log,prot_align,'>',prot_tree]
            os.system(' '.join(cmd))
        else:
            cmd = [FASTTREE,'-seed','12345','-boot',str(args.boot),'-gamma','-log',fasttree_log,prot_align,
            '>',prot_tree]
            os.system(' '.join(cmd))
    if args.type == 'nuc':
        if args.fasttree_args:
            fast_args = ''.join(args.fasttree_args).split(',')
            fast_args = [a.replace(a,'-'+a) for a in fast_args]
            cmd = [FASTTREE,'-nt','-seed','12345','-gamma','-boot',str(args.boot),' '.join(fast_args),'-log',
            fasttree_log,nuc_align,'>',nuc_tree]
            os.system(' '.join(cmd))
        else:
            cmd = [FASTTREE,'-nt','-seed','12345','-boot',str(args.boot),'-gamma','-log',fasttree_log,nuc_align,
            '>',nuc_tree]
            os.system(' '.join(cmd))
    if args.type == 'both':
        if args.fasttree_args:
            fast_args = ''.join(args.fasttree_args).split(',')
            fast_args = [a.replace(a,'-'+a) for a in fast_args]
            cmd = [FASTTREE,'-seed','12345','-gamma','-boot',str(args.boot),' '.join(fast_args),'-log',
            fasttree_log,prot_align,'>',prot_tree]
            os.system(' '.join(cmd))
            cmd = [FASTTREE,'-nt','-seed','12345','-gamma','-boot',str(args.boot),' '.join(fast_args),'-log',
            fasttree_log,nuc_align,'>',nuc_tree]
            os.system(' '.join(cmd))
        else:
            cmd = [FASTTREE,'-seed','12345','-boot',str(args.boot),'-gamma','-log',fasttree_log,prot_align,
            '>',prot_tree]
            os.system(' '.join(cmd))
            cmd = [FASTTREE,'-nt','-seed','12345','-boot',str(args.boot),'-gamma','-log',fasttree_log,nuc_align,
            '>',nuc_tree]
            os.system(' '.join(cmd))
if args.phylogeny == 'raxml':
    if args.type == 'prot':
        if not args.raxml_threads:
            if args.raxml_outgroup:
                if args.raxml_args:
                    rax_args = ''.join(args.raxml_args).split(',')
                    rax_args = [a.replace(a,'-'+a) for a in rax_args]
                    cmd = [RAXML+'-PTHREADS','-o',''.join(args.raxml_outgroup),' '.join(rax_args),
                    '-m',args.raxml_prot_sub,'-s',prot_align,'-x','12345','-p','12345','-#',str(args.boot),
                    '-k','-T',str(args.cpus),'-f',args.raxml_algorithm,'-n','prot_out']
                    os.system(' '.join(cmd))
                else:
                    cmd = [RAXML+'-PTHREADS','-o',''.join(args.raxml_outgroup),'-m',args.raxml_prot_sub,'-s',
                    prot_align,'-x','12345','-p','12345','-#',str(args.boot),'-k','-T',str(args.cpus),'-f',
                    args.raxml_algorithm,'-n','prot_out']
                    os.system(' '.join(cmd))
            elif args.raxml_args:
                rax_args = ''.join(args.raxml_args).split(',')
                rax_args = [a.replace(a,'-'+a) for a in rax_args]
                cmd = [RAXML+'-PTHREADS',' '.join(rax_args),'-m',args.raxml_prot_sub,'-s',prot_align,'-x','12345',
                '-p','12345','-#',str(args.boot),'-k','-T',str(args.cpus),'-f',args.raxml_algorithm,'-n','prot_out']
                os.system(' '.join(cmd))
            else:
                cmd = [RAXML+'-PTHREADS','-m',args.raxml_prot_sub,'-s',prot_align,'-x','12345','-p','12345',
                '-#',str(args.boot),'-k','-T',str(args.cpus),'-f',args.raxml_algorithm,'-n','prot_out']
                os.system(' '.join(cmd))
        if args.raxml_threads:
            if args.raxml_outgroup:
                if args.raxml_args:
                    rax_args = ''.join(args.raxml_args).split(',')
                    rax_args = [a.replace(a,'-'+a) for a in rax_args]
                    cmd = [RAXML,'-o',''.join(args.raxml_outgroup),' '.join(rax_args),
                    '-m',args.raxml_prot_sub,'-s',prot_align,'-x','12345','-p','12345','-#',str(args.boot),
                    '-k','-T',str(args.cpus),'-f',args.raxml_algorithm,'-n','prot_out']
                    os.system(' '.join(cmd))
                else:
                    cmd = [RAXML,'-o',''.join(args.raxml_outgroup),'-m',args.raxml_prot_sub,'-s',
                    prot_align,'-x','12345','-p','12345','-#',str(args.boot),'-k','-T',str(args.cpus),'-f',
                    args.raxml_algorithm,'-n','prot_out']
                    os.system(' '.join(cmd))
            elif args.raxml_args:
                rax_args = ''.join(args.raxml_args).split(',')
                rax_args = [a.replace(a,'-'+a) for a in rax_args]
                cmd = [RAXML,' '.join(rax_args),'-m',args.raxml_prot_sub,'-s',prot_align,'-x','12345',
                '-p','12345','-#',str(args.boot),'-k','-T',str(args.cpus),'-f',args.raxml_algorithm,'-n','prot_out']
                os.system(' '.join(cmd))
            else:
                cmd = [RAXML,'-m',args.raxml_prot_sub,'-s',prot_align,'-x','12345','-p','12345',
                '-#',str(args.boot),'-k','-T',str(args.cpus),'-f',args.raxml_algorithm,'-n','prot_out']
                os.system(' '.join(cmd))
        for file in os.listdir(currentdir):
            if file.startswith('RAxML_') and 'prot_out' in file:
                if 'bestTree' in file:
                    shutil.move(file,phylo_dir+'/{}_{}_{}_prot_tree.nwk'.format(
                    args.phylogeny,args.aligner,args.input))
                else:
                    shutil.move(file,phylo_dir)
    if args.type == 'nuc':
        if not args.raxml_threads:
            if args.raxml_outgroup:
                if args.raxml_args:
                    rax_args = ''.join(args.raxml_args).split(',')
                    rax_args = [a.replace(a,'-'+a) for a in rax_args]
                    cmd = [RAXML+'-PTHREADS','-o',''.join(args.raxml_outgroup),' '.join(rax_args),
                    '-m',args.raxml_nuc_sub,'-s',nuc_align,'-x','12345','-p','12345','-#',str(args.boot),
                    '-k','-T',str(args.cpus),'-f',args.raxml_algorithm,'-n','nuc_out']
                    os.system(' '.join(cmd))
                else:
                    cmd = [RAXML+'-PTHREADS','-o',''.join(args.raxml_outgroup),'-m',args.raxml_nuc_sub,'-s',
                    nuc_align,'-x','12345','-p','12345','-#',str(args.boot),'-k','-T',str(args.cpus),'-f',
                    args.raxml_algorithm,'-n','nuc_out']
                    os.system(' '.join(cmd))
            elif args.raxml_args:
                rax_args = ''.join(args.raxml_args).split(',')
                rax_args = [a.replace(a,'-'+a) for a in rax_args]
                cmd = [RAXML+'-PTHREADS',' '.join(rax_args),'-m',args.raxml_nuc_sub,'-s',nuc_align,'-x','12345',
                '-p','12345','-#',str(args.boot),'-k','-T',str(args.cpus),'-f',args.raxml_algorithm,'-n','nuc_out']
                os.system(' '.join(cmd))
            else:
                cmd = [RAXML+'-PTHREADS','-m',args.raxml_nuc_sub,'-s',nuc_align,'-x','12345','-p','12345',
                '-#',str(args.boot),'-k','-T',str(args.cpus),'-f',args.raxml_algorithm,'-n','nuc_out']
                os.system(' '.join(cmd))
        if args.raxml_threads:
            if args.raxml_outgroup:
                if args.raxml_args:
                    rax_args = ''.join(args.raxml_args).split(',')
                    rax_args = [a.replace(a,'-'+a) for a in rax_args]
                    cmd = [RAXML,'-o',''.join(args.raxml_outgroup),' '.join(rax_args),
                    '-m',args.raxml_nuc_sub,'-s',nuc_align,'-x','12345','-p','12345','-#',str(args.boot),
                    '-k','-T',str(args.cpus),'-f',args.raxml_algorithm,'-n','nuc_out']
                    os.system(' '.join(cmd))
                else:
                    cmd = [RAXML,'-o',''.join(args.raxml_outgroup),'-m',args.raxml_nuc_sub,'-s',
                    nuc_align,'-x','12345','-p','12345','-#',str(args.boot),'-k','-T',str(args.cpus),'-f',
                    args.raxml_algorithm,'-n','nuc_out']
                    os.system(' '.join(cmd))
            elif args.raxml_args:
                rax_args = ''.join(args.raxml_args).split(',')
                rax_args = [a.replace(a,'-'+a) for a in rax_args]
                cmd = [RAXML,' '.join(rax_args),'-m',args.raxml_nuc_sub,'-s',nuc_align,'-x','12345',
                '-p','12345','-#',str(args.boot),'-k','-T',str(args.cpus),'-f',args.raxml_algorithm,'-n','nuc_out']
                os.system(' '.join(cmd))
            else:
                cmd = [RAXML,'-m',args.raxml_nuc_sub,'-s',nuc_align,'-x','12345','-p','12345',
                '-#',str(args.boot),'-k','-T',str(args.cpus),'-f',args.raxml_algorithm,'-n','nuc_out']
                os.system(' '.join(cmd))
        for file in os.listdir(currentdir):
            if file.startswith('RAxML_') and 'nuc_out' in file:
                if 'bestTree' in file:
                    shutil.move(file,phylo_dir+'/{}_{}_{}_nuc_tree.nwk'.format(
                    args.phylogeny,args.aligner,args.input))
                else:
                    shutil.move(file,phylo_dir)
    if args.type == 'both':
        if not args.raxml_threads:
            if args.raxml_outgroup:
                if args.raxml_args:
                    rax_args = ''.join(args.raxml_args).split(',')
                    rax_args = [a.replace(a,'-'+a) for a in rax_args]
                    cmd = [RAXML+'-PTHREADS','-o',''.join(args.raxml_outgroup),' '.join(rax_args),
                    '-m',args.raxml_prot_sub,'-s',prot_align,'-x','12345','-p','12345','-#',str(args.boot),
                    '-k','-T',str(args.cpus),'-f',args.raxml_algorithm,'-n','prot_out']
                    os.system(' '.join(cmd))
                    rax_args = ''.join(args.raxml_args).split(',')
                    rax_args = [a.replace(a,'-'+a) for a in rax_args]
                    cmd = [RAXML+'-PTHREADS','-o',''.join(args.raxml_outgroup),' '.join(rax_args),
                    '-m',args.raxml_nuc_sub,'-s',nuc_align,'-x','12345','-p','12345','-#',str(args.boot),
                    '-k','-T',str(args.cpus),'-f',args.raxml_algorithm,'-n','nuc_out']
                    os.system(' '.join(cmd))
                else:
                    cmd = [RAXML+'-PTHREADS','-o',''.join(args.raxml_outgroup),'-m',args.raxml_prot_sub,'-s',
                    prot_align,'-x','12345','-p','12345','-#',str(args.boot),'-k','-T',str(args.cpus),'-f',
                    args.raxml_algorithm,'-n','prot_out']
                    os.system(' '.join(cmd))
                    cmd = [RAXML+'-PTHREADS','-o',''.join(args.raxml_outgroup),'-m',args.raxml_nuc_sub,'-s',
                    nuc_align,'-x','12345','-p','12345','-#',str(args.boot),'-k','-T',str(args.cpus),'-f',
                    args.raxml_algorithm,'-n','nuc_out']
                    os.system(' '.join(cmd))
            elif args.raxml_args:
                rax_args = ''.join(args.raxml_args).split(',')
                rax_args = [a.replace(a,'-'+a) for a in rax_args]
                cmd = [RAXML+'-PTHREADS',' '.join(rax_args),'-m',args.raxml_prot_sub,'-s',prot_align,'-x','12345',
                '-p','12345','-#',str(args.boot),'-k','-T',str(args.cpus),'-f',args.raxml_algorithm,'-n','prot_out']
                os.system(' '.join(cmd))
                rax_args = ''.join(args.raxml_args).split(',')
                rax_args = [a.replace(a,'-'+a) for a in rax_args]
                cmd = [RAXML+'-PTHREADS',' '.join(rax_args),'-m',args.raxml_nuc_sub,'-s',nuc_align,'-x','12345',
                '-p','12345','-#',str(args.boot),'-k','-T',str(args.cpus),'-f',args.raxml_algorithm,'-n','nuc_out']
                os.system(' '.join(cmd))
            else:
                cmd = [RAXML+'-PTHREADS','-m',args.raxml_prot_sub,'-s',prot_align,'-x','12345','-p','12345',
                '-#',str(args.boot),'-k','-T',str(args.cpus),'-f',args.raxml_algorithm,'-n','prot_out']
                os.system(' '.join(cmd))
                cmd = [RAXML+'-PTHREADS','-m',args.raxml_nuc_sub,'-s',nuc_align,'-x','12345','-p','12345',
                '-#',str(args.boot),'-k','-T',str(args.cpus),'-f',args.raxml_algorithm,'-n','nuc_out']
                os.system(' '.join(cmd))
        if args.raxml_threads:
            if args.raxml_outgroup:
                if args.raxml_args:
                    rax_args = ''.join(args.raxml_args).split(',')
                    rax_args = [a.replace(a,'-'+a) for a in rax_args]
                    cmd = [RAXML,'-o',''.join(args.raxml_outgroup),' '.join(rax_args),
                    '-m',args.raxml_prot_sub,'-s',prot_align,'-x','12345','-p','12345','-#',str(args.boot),
                    '-k','-T',str(args.cpus),'-f',args.raxml_algorithm,'-n','prot_out']
                    os.system(' '.join(cmd))
                    rax_args = ''.join(args.raxml_args).split(',')
                    rax_args = [a.replace(a,'-'+a) for a in rax_args]
                    cmd = [RAXML,'-o',''.join(args.raxml_outgroup),' '.join(rax_args),
                    '-m',args.raxml_nuc_sub,'-s',nuc_align,'-x','12345','-p','12345','-#',str(args.boot),
                    '-k','-T',str(args.cpus),'-f',args.raxml_algorithm,'-n','nuc_out']
                    os.system(' '.join(cmd))
                else:
                    cmd = [RAXML,'-o',''.join(args.raxml_outgroup),'-m',args.raxml_prot_sub,'-s',
                    prot_align,'-x','12345','-p','12345','-#',str(args.boot),'-k','-T',str(args.cpus),'-f',
                    args.raxml_algorithm,'-n','prot_out']
                    os.system(' '.join(cmd))
                    cmd = [RAXML,'-o',''.join(args.raxml_outgroup),'-m',args.raxml_nuc_sub,'-s',
                    nuc_align,'-x','12345','-p','12345','-#',str(args.boot),'-k','-T',str(args.cpus),'-f',
                    args.raxml_algorithm,'-n','nuc_out']
                    os.system(' '.join(cmd))
            elif args.raxml_args:
                rax_args = ''.join(args.raxml_args).split(',')
                rax_args = [a.replace(a,'-'+a) for a in rax_args]
                cmd = [RAXML,' '.join(rax_args),'-m',args.raxml_prot_sub,'-s',prot_align,'-x','12345',
                '-p','12345','-#',str(args.boot),'-k','-T',str(args.cpus),'-f',args.raxml_algorithm,'-n','prot_out']
                os.system(' '.join(cmd))
                rax_args = ''.join(args.raxml_args).split(',')
                rax_args = [a.replace(a,'-'+a) for a in rax_args]
                cmd = [RAXML,' '.join(rax_args),'-m',args.raxml_nuc_sub,'-s',nuc_align,'-x','12345',
                '-p','12345','-#',str(args.boot),'-k','-T',str(args.cpus),'-f',args.raxml_algorithm,'-n','nuc_out']
                os.system(' '.join(cmd))
            else:
                cmd = [RAXML,'-m',args.raxml_prot_sub,'-s',prot_align,'-x','12345','-p','12345',
                '-#',str(args.boot),'-k','-T',str(args.cpus),'-f',args.raxml_algorithm,'-n','prot_out']
                os.system(' '.join(cmd))
                cmd = [RAXML,'-m',args.raxml_nuc_sub,'-s',nuc_align,'-x','12345','-p','12345',
                '-#',str(args.boot),'-k','-T',str(args.cpus),'-f',args.raxml_algorithm,'-n','nuc_out']
                os.system(' '.join(cmd))
        for file in os.listdir(currentdir):
            if file.startswith('RAxML_') and 'prot_out' in file:
                if 'bestTree' in file:
                    shutil.move(file,phylo_dir+'/{}_{}_{}_prot_tree.nwk'.format(
                    args.phylogeny,args.aligner,args.input))
                else:
                    shutil.move(file,phylo_dir)
        for file in os.listdir(currentdir):
            if file.startswith('RAxML_') and 'nuc_out' in file:
                if 'bestTree' in file:
                    shutil.move(file,phylo_dir+'/{}_{}_{}_nuc_tree.nwk'.format(
                    args.phylogeny,args.aligner,args.input))
                else:
                    shutil.move(file,phylo_dir)

print('done')
