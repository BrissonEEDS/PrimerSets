import click
import subprocess
from click._compat import filename_to_ui
import os, stat
from pyfaidx import Fasta
from swga.utils.resources import get_swga_opts
from swga.utils.options import cfg_from_opts

welcome_message = """
## SWGA Initialization ---------------------------

This will set up a workspace for SWGA in the current directory ({cwd}).

As part of this process, a default parameters file will be created that contains
some reasonable options for each part of the pipeline. This will be
named "{default_parameters_name}" and should be modified as appropriate for your
particular project. Please note that you should only with edit this file with a
plain-text editor like Notepad, TextEdit, or gedit; or with command-line tools
such as nano, vim, or emacs. 

The values specified in {default_parameters_name} will be used if the
corresponding values are not specified on the command line.
"""

fg_message = """
Foreground genome: {fg_genome_fp}
  Length:  {fg_length} bp
  Records: {fg_nrecords}
"""
bg_message = """
Background genome: {bg_genome_fp}
  Length:  {bg_length} bp
  Records: {bg_nrecords}
"""

finish_message = """Done!"""


@click.command()
@click.option("-f", "--fg_genome_fp",
              type=click.Path(exists=True, resolve_path=True))
@click.option("-b", "--bg_genome_fp",
              type=click.Path(exists=True, resolve_path=True))
@click.option("-e", "--exclude_fp",
              type=click.Path(exists=True, resolve_path=True))
def main(fg_genome_fp, bg_genome_fp, exclude_fp):

    default_parameters_name = "parameters.cfg"
    cwd = os.getcwd()

    # Print the welcome message ------------------------------
    click.echo(click.style(welcome_message.format(**locals()),
                           fg = "blue"))

    # Prompt for the foreground genome, if not already specified                                   
    if (not fg_genome_fp):
        fg_genome_fp = click.prompt("Enter path to foreground genome file, in "
                                    "FASTA format", 
                                    type=click.Path(exists=True,
                                                    resolve_path=True))
    fg_length, fg_nrecords = fasta_stats(fg_genome_fp)
    click.echo(click.style(fg_message.format(**locals()), fg="green"))


    # Prompt for background genome, if not already specified
    if (not bg_genome_fp):
        bg_genome_fp = click.prompt("Enter path to background genome file, in "
                                    "FASTA format", 
                                    type=click.Path(exists=True, 
                                                    resolve_path=True))
    bg_length, bg_nrecords = fasta_stats(bg_genome_fp)
    click.echo(click.style(bg_message.format(**locals()), fg="green"))


    # Prompt for a file containing sequences to exclude, if not already given
    if (not exclude_fp):
        if click.confirm(
                "Do you want to add a FASTA file that will be used to exclude "
                "primers? For instance, to avoid primers binding to the "
                "mitochondrial genome, add the path to that sequence here."):
            
            exclude_fp = click.prompt(
                "Enter path to exclusionary sequence(s), in FASTA format",
                type=click.Path(exists=True, resolve_path=True))

            click.echo(click.style(
                "Exclusionary sequences file: {}".format(exclude_fp), fg="red"))
        else:
            click.echo(
                click.style("No exclusionary sequences specified.", fg="green"))
            exclude_fp = ""


    opts = get_swga_opts()
    default_parameters = cfg_from_opts(opts)

    cfg_fp = os.path.join(cwd, default_parameters_name)
    if os.path.isfile(cfg_fp):
        click.confirm("Existing file `%s` will be overwritten. Continue?" 
                      % default_parameters_name, abort=True)

    # Calculate default binding frequencies for foreground and background
    # Fg is based on a binding rate of 1/100000 bp/binding site
    # Fg is based on a binding rate of 1/150000 bp/binding site
    min_fg_rate = 0.00001
    min_fg_bind = int(min_fg_rate*float(fg_length))
    max_bg_rate = 0.0000067
    max_bg_bind = int(max_bg_rate*float(bg_length))
    with open(os.path.join(cwd, default_parameters_name), "wb") as cfg_file:
        cfg_file.write(default_parameters.format(**locals()))
    
    click.echo(finish_message)
              

def fasta_stats(fasta_fp):
    check_empty_lines(fasta_fp)
    try:
        fasta = Fasta(fasta_fp)
        length = sum([len(r) for r in fasta])
        nrecords = len(fasta.keys())
        return length, nrecords
    except:
        click.echo(click.style("\nError reading %s: invalid FASTA format?" %
                               fasta_fp, fg = "red"))
        raise


def check_empty_lines(fasta_fp):
    '''
    The pyfaidx module has issues with blank lines in Fasta files. This checks
    for their presence and offers to remove them in-place.
    '''
    basename = os.path.basename(fasta_fp)
    try:
        check = subprocess.check_output(
            "grep -n '^$' {}".format(fasta_fp), shell=True)
    except subprocess.CalledProcessError:
        check = None
    if check:
        click.confirm(
            "`{}` has blank lines, which can interfere with SWGA."
            " Is it okay to remove these lines from the file?"
            .format(basename), abort=True)
        # sed -ie gets around the fact that BSD sed -i requires a backup extension
        # with sed, but we don't want to create a backup. Ensures linux/osx
        # compatibility. 
        subprocess.check_call("sed -ie '/^$/d' \"{}\"".format(fasta_fp),
                              shell=True) 
            

def monkeypatch_method(cls):
    def decorator(func):
        setattr(cls, func.__name__, func)
        return func
    return decorator

@monkeypatch_method(click.Path)
def convert(self, value, param, ctx):
    value = value.strip()
    rv = value
    if self.resolve_path:
        rv = os.path.realpath(rv)

    try:
        st = os.stat(rv)
    except OSError:
        if not self.exists:
            return rv
        self.fail('%s "%s" does not exist.' % (
            self.path_type,
            filename_to_ui(value)
        ), param, ctx)

    if not self.file_okay and stat.S_ISREG(st.st_mode):
        self.fail('%s "%s" is a file.' % (
            self.path_type,
            filename_to_ui(value)
        ), param, ctx)
    if not self.dir_okay and stat.S_ISDIR(st.st_mode):
        self.fail('%s "%s" is a directory.' % (
            self.path_type,
            filename_to_ui(value)
        ), param, ctx)
    if self.writable and not os.access(value, os.W_OK):
        self.fail('%s "%s" is not writable.' % (
            self.path_type,
            filename_to_ui(value)
        ), param, ctx)
    if self.readable and not os.access(value, os.R_OK):
        self.fail('%s "%s" is not readable.' % (
            self.path_type,
            filename_to_ui(value)
        ), param, ctx)

    return rv

    


if __name__ == "__main__":
    main()
