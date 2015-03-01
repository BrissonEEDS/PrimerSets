# -*- coding: utf-8 -*-
'''
Selects valid primers from a database according to various criteria.

This code could be a lot more succinct, but it's largely broken apart
to provide useful messages during operation (like how many primers satisfy
each criteria) and for performance reasons (i.e. we only calculate melt
temps on primers that pass the binding rate thresholds, and only calculate
binding locations on primers that pass all filters.)

Since all the results of the calculations are stored in the primer db, 
the user can run this command multiple times to tune parameters and
subsequent runs will be much faster.
'''
import json

import swga.locate as locate
import swga.primers
from swga.clint.textui import progress
from swga.commands import Command
from swga.commands.activate import update_tms, update_locations
from swga.database import Primer, update_in_chunks
from swga.melting import Tm


def main(argv, cfg_file):
    cmd = Command('filter', cfg_file=cfg_file)
    cmd.parse_args(argv)

    swga.database.init_db(cmd.primer_db)
    
    # If we have an input file, use that. Otherwise pull from db
    if cmd.input:
        with open(cmd.input, 'rb') as infile:
            primers = swga.primers.read_primer_list(
                infile,
                cmd.fg_genome_fp,
                cmd.bg_genome_fp)  
    else:
        cmd.skip_filtering = False
        primers = Primer.select()

    # Undo all active marks, if any
    deactivate_all_primers()
    
    if not cmd.skip_filtering:
        primers = filter_primers(
            primers,
            cmd.min_fg_bind,
            cmd.max_bg_bind,
            cmd.fg_length,
            cmd.bg_length,
            cmd.min_tm,
            cmd.max_tm,
            cmd.max_primers)
    
    update_locations(primers, cmd.fg_genome_fp)
    n_active = activate_primers(primers)
    if int(n_active) < int(cmd.max_primers):
    	swga.warn(
    		"Fewer than {} primers were selected. Only {} passed all "
		"the filters You may want to try less "
    		"restrictive filtering parameters.".format(cmd.max_primers, n_active))


def deactivate_all_primers():
    """Resets all active marks on primers."""
    Primer.update(active=False).execute()

    
def activate_primers(primers):
    """
    Marks as active all the kmers passed to it.
    """
    n_active = Primer.update(active=True).where(Primer.seq << primers).execute()
    swga.message("Marked {} primers as active.".format(n_active))
    return n_active


def filter_primers(
        primers,
        min_fg_bind,
        max_bg_bind,
        fg_length,
        bg_length,
        min_tm,
        max_tm,
        max_primers):
    """
    Takes a list of sequences and retrieves them in the database, then returns
    those sequences that pass various criteria.
    """
    primers = Primer.select().where(Primer.seq << primers)
    fg_min_freq = min_fg_bind
    bg_max_freq = max_bg_bind

    # Find primers that pass the binding rate thresholds
    fgp = Primer.select().where((Primer.seq << primers) &
                                (Primer.fg_freq >= fg_min_freq))
    swga.message("{} primers bind foreground genome with freq >= {} sites"
                 .format(fgp.count(), min_fg_bind))
    
    bgp = Primer.select().where((Primer.seq << primers) &
                                (Primer.bg_freq <= bg_max_freq))
    swga.message("{} primers bind background genome with freq <= {} sites"
                 .format(bgp.count(), max_bg_bind))

    candidates = Primer.select().where((Primer.seq << primers) &
                                       (Primer.seq << fgp) &
                                       (Primer.seq << bgp))
    swga.message("{} primers pass both fg and bg binding freq filters"
                .format(candidates.count()))

    # Add melt temp for any primer that doesn't have it yet
    update_tms(candidates)
    
    valid_primers = Primer.select().where((Primer.seq << candidates) &
                                          (Primer.tm <= max_tm) &
                                          (Primer.tm >= min_tm))
    swga.message("{} of those primers have a melting temp within given range"
                 .format(valid_primers.count()))

    # Sort by background binding rate (smallest -> largest) and select top `n`,
    # then sort those by ratio (highest -> lowest)
    first_pass = (Primer.select()
                  .where(Primer.seq << valid_primers)
                  .order_by(Primer.bg_freq)
                  .limit(max_primers))

    second_pass = (Primer.select()
                   .where(Primer.seq << first_pass)
                   .order_by(Primer.ratio.desc()))

    return second_pass
