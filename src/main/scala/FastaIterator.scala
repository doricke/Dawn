
import scala.collection.mutable.{ArrayBuffer, Map}

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

// ****************************************************************************
class FastaIterator( fileName: String ) extends InputFile( fileName ) {
  
  private[this] var fasta: FastaSequence = new FastaSequence()

  // **************************************************************************
  def getFasta: FastaSequence = fasta

  // **************************************************************************
  def fastas(): ArrayBuffer[FastaSequence] = {
    val seqs: ArrayBuffer[FastaSequence] = ArrayBuffer()

    // Read in all of the sequences
    while ( endOfFile == false ) {
      nextFasta( "" )
      seqs += fasta
    }  // while

    seqs
  }  // fastas

  // **************************************************************************
  def fastasToMap: Map[String, FastaSequence] = {
    val seqs = collection.mutable.Map[String, FastaSequence]()

    // Read in all of the sequences
    while ( endOfFile == false ) {
      nextFasta( "" )

      if ( fasta.sequence.length > 0 )
        seqs += (fasta.name -> fasta)
    }  // while

    seqs
  }  // fastasToMap

  // **************************************************************************
  def seqsToMap( species: String, control: String ): Map[String, String] = {
    val seqs = collection.mutable.Map[String, String]()

    // Read in all of the sequences
    while ( endOfFile == false ) {
      nextFasta( "" )
      if ( ( fasta.description contains species ) || ( fasta.description contains control ) ) {
        fasta.sequence = fasta.sequence.replace( 'U', 'T' )
        seqs += (fasta.name -> fasta.sequence)
        // println( "miRNA: " + fasta.name )
      }  // if
    }  // while

    seqs
  }  // seqsToMap

  // **************************************************************************
  def nextFasta( barcode: String ): FastaSequence = {
    fasta = new FastaSequence()

    // println( "FastaIterator.next: line = " + line )

    fasta.parseHeader( line )

    readSequence()
    fasta.mid = barcode
    fasta
  }  // method next

  // **************************************************************************
  def countSequences: Int = {
    var count: Int = 0

    nextLine
    while ( endOfFile == false ) {
      if ( line.charAt( 0 ) == '>' ) 
        count += 1
      nextLine
    }  // while

    count
  }  // method countSequences


  // **************************************************************************
  def readSequence() {
    var seq: String = ""

    nextLine
    while ( ( endOfFile == false ) && ( line.length > 0 ) && ( line.charAt( 0 ) != '>' ) )
    {
      if ( line.charAt( 0 ) != '>' )
      {
        seq += line
        nextLine
      }  // if
    }  // while

    // println( "FastaIterator.readSequence: seq = " + seq )

    fasta.sequence = seq
  }  // method readSequence


  // **************************************************************************

}  // class FastaIterator
