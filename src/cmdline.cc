#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef FIX_UNUSED
#define FIX_UNUSED(X) (void) (X) /* avoid warnings for unused params */
#endif

#include <getopt.h>

#include "cmdline.hh"

std::string input_struct;
std::string input_file;
std::string output_file;
int subopt;

static char *package_name = 0;

const char *args_info_purpose = "Possibly pseudoknotted minimum free energy folding of RNAs";

const char *args_info_usage = "Usage: Iterative-HFold"  "[options] [sequence]";

const char *args_info_versiontext = "";

const char *args_info_description = "Read RNA sequence from stdin or cmdline; predict minimum\nfree energy and optimum structure using the possibly pseuodknotted\nMFE RNA folding algorithm of Jabbari and Condon, 2014.";

const char *args_info_help[] = {
  "  -h, --help             Print help and exit",
  "  -V, --version          Print version and exit",
  "  -v, --verbose          Give verbose output",
  "  -r, --input-structure  Give a restricted structure as an input structure",
  "  -i, --input-file       Give a path to an input file containing the sequence (and input structure if known)",
  "  -o  --output-file      Give a path to an output file which will the sequence, and its structure and energy",
  "  -n, --opt              Specify the number of suboptimal structures to output (default is 1)",
  "  -p  --pk-free          Specify whether you only want the pseudoknot-free simfold structure to be calculated"

  "\nThe input sequence is read from standard input, unless it is\ngiven on the command line.\n",
  
};

static void clear_given (struct args_info *args_info);
static void clear_args (struct args_info *args_info);

static int cmdline_parser_internal (int argc, char **argv, struct args_info *args_info, const char *additional_error);

typedef enum {ARG_NO} cmdline_parser_arg_type;

static char *gengetopt_strdup (const char *s);

static void init_args_info(struct args_info *args_info)
{
  args_info->help_help = args_info_help[0] ;
  args_info->version_help = args_info_help[1] ;
  args_info->verbose_help = args_info_help[2] ;
  args_info->input_structure_help = args_info_help[3] ;
  args_info->input_file_help = args_info_help[4] ;
  args_info->output_file_help = args_info_help[5] ;
  args_info->subopt_help = args_info_help[7] ;
  args_info->pk_free_help = args_info_help[8] ;


  
}
void
cmdline_parser_print_version (void)
{
  printf (" %s\n",(strlen(package_name) ? package_name : "Iterative-HFold"));

  printf ("%s %s\n",
     (strlen(Iterative_HFold_CMDLINE_PACKAGE_NAME) ? Iterative_HFold_CMDLINE_PACKAGE_NAME : "Iterative-HFold"),
     Iterative_HFold_CMDLINE_VERSION);

  if (strlen(args_info_versiontext) > 0)
    printf("\n%s\n", args_info_versiontext);
}

static void print_help_common(void)
{
	size_t len_purpose = strlen(args_info_purpose);
	size_t len_usage = strlen(args_info_usage);

	if (len_usage > 0) {
		printf("%s\n", args_info_usage);
	}
	if (len_purpose > 0) {
		printf("%s\n", args_info_purpose);
	}

	if (len_usage || len_purpose) {
		printf("\n");
	}

	if (strlen(args_info_description) > 0) {
		printf("%s\n\n", args_info_description);
	}
}
void cmdline_parser_print_help (void){
  
  print_help_common();
  int i = 0;
  int end = sizeof(args_info_help)/sizeof(args_info_help[0]);
  while (i<end) printf("%s\n", args_info_help[i++]);
}

static void clear_given (struct args_info *args_info)
{
  args_info->help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->verbose_given = 0 ;
  args_info->input_structure_given = 0 ;
  args_info->input_file_given = 0 ;
  args_info->output_file_given = 0 ;
  args_info->subopt_given = 0 ;
  args_info->pk_free_given = 0 ;
}

static void clear_args (struct args_info *args_info)
{
  FIX_UNUSED (args_info);
  
}

static void cmdline_parser_release (struct args_info *args_info)
{
  unsigned int i;
  
  
  for (i = 0; i < args_info->inputs_num; ++i)
    free (args_info->inputs [i]);

  if (args_info->inputs_num)
    free (args_info->inputs);

  clear_given (args_info);
}

void cmdline_parser_init (struct args_info *args_info)
{
  clear_given (args_info);
  clear_args (args_info);
  init_args_info (args_info);

  args_info->inputs = 0;
  args_info->inputs_num = 0;
}

void cmdline_parser_free (struct args_info *args_info)
{
  cmdline_parser_release (args_info);
}

/** @brief replacement of strdup, which is not standard */
char *
gengetopt_strdup (const char *s)
{
  char *result = 0;
  if (!s)
    return result;

  result = (char*)malloc(strlen(s) + 1);
  if (result == (char*)0)
    return (char*)0;
  strcpy(result, s);
  return result;
}




int cmdline_parser (int argc, char **argv, struct args_info *args_info)
{
  return cmdline_parser2 (argc, argv, args_info);
}

int cmdline_parser2 (int argc, char **argv, struct args_info *args_info)
{
  int result;

  result = cmdline_parser_internal (argc, argv, args_info, 0);

  if (result == EXIT_FAILURE)
    {
      cmdline_parser_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

/**
 * @brief updates an option
 * @param field the generic pointer to the field to update
 * @param orig_field the pointer to the orig field
 * @param field_given the pointer to the number of occurrence of this option
 * @param prev_given the pointer to the number of occurrence already seen
 * @param value the argument for this option (if null no arg was specified)
 * @param possible_values the possible values for this option (if specified)
 * @param default_value the default value (in case the option only accepts fixed values)
 * @param arg_type the type of this option
 * @param check_ambiguity @see cmdline_parser_params.check_ambiguity
 * @param override @see cmdline_parser_params.override
 * @param no_free whether to free a possible previous value
 * @param multiple_option whether this is a multiple option
 * @param long_opt the corresponding long option
 * @param short_opt the corresponding short option (or '-' if none)
 * @param additional_error possible further error specification
 */
static
int update_arg(void *field, char **orig_field,unsigned int *field_given,
               unsigned int *prev_given, char *value, const char *possible_values[],
               const char *default_value,cmdline_parser_arg_type arg_type,
               int no_free, int multiple_option,
               const char *long_opt, char short_opt,const char *additional_error)
{
  char *stop_char = 0;
  const char *val = value;
  int found;
  FIX_UNUSED (field);

  stop_char = 0;
  found = 0;

  if (!multiple_option && prev_given && (*prev_given || (0 && *field_given)))
    {
      if (short_opt != '-')
        fprintf (stderr, "%s: `--%s' (`-%c') option given more than once%s\n", 
               package_name, long_opt, short_opt,
               (additional_error ? additional_error : ""));
      else
        fprintf (stderr, "%s: `--%s' option given more than once%s\n", 
               package_name, long_opt,
               (additional_error ? additional_error : ""));
      return 1; /* failure */
    }

  FIX_UNUSED (default_value);
    
  if (field_given && *field_given && ! 0)
    return 0;
  if (prev_given)
    (*prev_given)++;
  if (field_given)
    (*field_given)++;
  if (possible_values)
    val = possible_values[found];

  switch(arg_type) {
  default:
    break;
  };

	FIX_UNUSED(stop_char);
	FIX_UNUSED(val);
	
  /* store the original value */
  switch(arg_type) {
  case ARG_NO:
    break;
  default:
    if (value && orig_field) {
      if (no_free) {
        *orig_field = value;
      } else {
        if (*orig_field)
          free (*orig_field); /* free previous string */
        *orig_field = gengetopt_strdup (value);
      }
    }
  };

  return 0; /* OK */
}



int cmdline_parser_internal (int argc, char **argv, struct args_info *args_info, const char *additional_error){
  int c;	/* Character of the parsed option.  */

  int error_occurred = 0;
  struct args_info local_args_info;
  
 
  
  package_name = argv[0];


  cmdline_parser_init (args_info);

  cmdline_parser_init (&local_args_info);

  optarg = 0;
  optind = 0;
  optopt = '?';

  while (1)
    {
      int option_index = 0;

      static struct option long_options[] = {
        { "help",	0, NULL, 'h' },
        { "version",	0, NULL, 'V' },
        { "verbose",	0, NULL, 'v' },
        { "input-structure",	required_argument, NULL, 'r' },
        { "input-file",	required_argument, NULL, 'i' },
        { "ouput-file",	required_argument, NULL, 'o' },
        { "subopt",	required_argument, NULL, 'n' },
        { "pk-free",	0, NULL, 'p' },
        { 0,  0, 0, 0 }
      };

      c = getopt_long (argc, argv, "hVvr:i:o:n:p", long_options, &option_index);

      if (c == -1) break;	/* Exit from `while (1)' loop.  */

      switch (c)
        {
        case 'h':	/* Print help and exit.  */
          cmdline_parser_print_help ();
          cmdline_parser_free (&local_args_info);
          exit (EXIT_SUCCESS);

        case 'V':	/* Print version and exit.  */
          cmdline_parser_print_version ();
          cmdline_parser_free (&local_args_info);
          exit (EXIT_SUCCESS);
        
        case 'v':	/* Turn on verbose output.  */
        
        
          if (update_arg( 0 , 
               0 , &(args_info->verbose_given),
              &(local_args_info.verbose_given), optarg, 0, 0, ARG_NO,0, 0,"verbose", 'v',additional_error))
            goto failure;
        
          break;

          case 'n':	/* Specify number of outputs.  */
        
        
          if (update_arg( 0 , 
               0 , &(args_info->subopt_given),
              &(local_args_info.subopt_given), optarg, 0, 0, ARG_NO,0, 0,"subopt", 'n',additional_error))
            goto failure;

            subopt = strtol(optarg,NULL,10);
        
          break;
        case 'r':	/* Specify restricted structure.  */
        
        
          if (update_arg( 0 , 
               0 , &(args_info->input_structure_given),
              &(local_args_info.input_structure_given), optarg, 0, 0, ARG_NO,0, 0,"input-structure", 'r',additional_error))
            goto failure;

            input_struct = optarg;
        
          break;
          case 'i':	/* Specify input file.  */
        
        
          if (update_arg( 0 , 
               0 , &(args_info->input_file_given),
              &(local_args_info.input_file_given), optarg, 0, 0, ARG_NO,0, 0,"input-file", 'i',additional_error))
            goto failure;

            input_file = optarg;
        
          break;
          case 'o':	/* Specify output file.  */
        
        
          if (update_arg( 0 , 
               0 , &(args_info->output_file_given),
              &(local_args_info.output_file_given), optarg, 0, 0, ARG_NO,0, 0,"output-file", 'o',additional_error))
            goto failure;

            output_file = optarg;
        
          break;

          case 'p':	/* Specify if pseudoknot-free.  */
        
        
          if (update_arg( 0 , 
               0 , &(args_info->pk_free_given),
              &(local_args_info.pk_free_given), optarg, 0, 0, ARG_NO,0, 0,"pk-free", 'p',additional_error))
            goto failure;
        
          break;

        case 0:	/* Long option with no short option */
          
          break;
        case '?':	/* Invalid option.  */
          /* `getopt_long' already printed an error message.  */
          goto failure;

        default:	/* bug: option not considered.  */
          fprintf (stderr, "%s: option unknown: %c%s\n", package_name, c, (additional_error ? additional_error : ""));
          abort ();
        } /* switch */
    } /* while */


  cmdline_parser_release (&local_args_info);

  if ( error_occurred )
    return (EXIT_FAILURE);

  if (optind < argc)
    {
      int i = 0 ;
      int found_prog_name = 0;
      /* whether program name, i.e., argv[0], is in the remaining args
         (this may happen with some implementations of getopt,
          but surely not with the one included by gengetopt) */

      i = optind;
      while (i < argc)
        if (argv[i++] == argv[0]) {
          found_prog_name = 1;
          break;
        }
      i = 0;

      args_info->inputs_num = argc - optind - found_prog_name;
      args_info->inputs =
        (char **)(malloc ((args_info->inputs_num)*sizeof(char *))) ;
      while (optind < argc)
        if (argv[optind++] != argv[0])
          args_info->inputs[ i++ ] = gengetopt_strdup (argv[optind-1]) ;
    }

  return 0;

  failure:
  
  cmdline_parser_release (&local_args_info);
  return (EXIT_FAILURE);
}