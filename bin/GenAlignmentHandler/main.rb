#!/usr/bin/ruby
#
#
# Takes the genthreader alignments and outputs jalview alignment and annotation
# files and a contact residue consensus
#

require "./bin/GenAlignmentHandler/lib/gen_alignment_handler.rb"

name = ARGV[0]
csa = ARGV[1]

gen_alignment_handler = GenAlignmentHandler.new(name,csa)
gen_alignment_handler.readPresults
gen_alignment_handler.readAlignments
gen_alignment_handler.happen
