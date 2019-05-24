
import scala.collection.mutable.{ArrayBuffer, Map, Set}

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
class Protein(val fasta: FastaSequence, val wordSize: Int) {
  val seqTools = new SeqTools()

  val mers: collection.mutable.Set[String] = seqTools.seqSplit( fasta.sequence, wordSize )

  val positions: collection.mutable.Map[String, ArrayBuffer[Int]] = seqTools.seqWords( fasta.sequence, wordSize )

  // **************************************************************************
}  // class Protein
