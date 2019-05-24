
import scala.collection.mutable.Map

//  Author:     Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
//  Copyright:  Copyright (c) 2015 Massachusetts Institute of Technology 
//  License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

class MsaSetup( queries: String, targets: String, pdb_names: String ) {

  // ***************************************************************************
  def longestSeq( seqs: Map[String, FastaSequence] ): FastaSequence =
  {
    // seqs.reduceLeft( _ max _._2.sequence.size )
    var longest: FastaSequence = seqs.head._2
    var length = longest.sequence.size

    seqs foreach {case (name, fasta) =>
      if ( fasta.sequence.length > length ) {
        longest = fasta
        length = fasta.sequence.length
      }  // if
    }  // foreach

    longest
  }  // longestSeq

  // ***************************************************************************
  def process(): Unit = {
    val wordSize = 3

    // Read in the targets.
    val targets_reader = if ( targets.size > 0 ) new FastaIterator( targets ) else new FastaIterator( queries )
    val target_seqs = targets_reader.fastasToMap

    // Map the targets to k-mers
    val target_mers = Map[String, Protein]()
    target_seqs foreach {case (name, fasta) =>
      target_mers += (name -> new Protein( fasta, wordSize ))
    }  // foreach

    val readers = Map[String, PdbReader]()

    // Set up the PDB structures.
    if ( pdb_names.size > 0 ) {
      // Read in the PDB names.
      val pdb_list = new InputFile( pdb_names )
      val pdbs = pdb_list.readContents().split( "\n" )

      pdbs foreach {case(pdb) =>
        val pdb_reader = new PdbReader( pdb )
        pdb_reader.parsePdb()
        val tokens = pdb.split( '.' )
        pdb_reader.chains foreach {case(chain, seq) =>
          if ( seq.length >= 20 ) {
            val fasta = new FastaSequence() 
            fasta.sequence = seq.toString
            fasta.name = tokens( 0 ) + "_" + chain
            target_mers += (fasta.name -> new Protein( fasta, wordSize ))
            readers += fasta.name -> pdb_reader
          }  // if
        }  // foreach
      }  // foreach
    }  // if

    // Set up the aligners.
    val aligner = new Aligner( queries, target_mers, wordSize, readers )

    // Distribute the queries to the aligners.
    if ( targets.size > 0 ) {
      val query_reader = new FastaIterator( queries )
      while ( query_reader.endOfFile == false ) {
        val fasta = query_reader.nextFasta( "" )
        // println( "Sending: " + fasta.name )
        aligner.process( fasta )
      }  // while
    }  else {
      aligner.process( longestSeq( target_seqs ) )
    }  // if

  }  // process

  // ***************************************************************************
}  // class MsaSetup
