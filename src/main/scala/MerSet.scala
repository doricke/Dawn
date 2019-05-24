
import scala.collection.mutable.{Map, ArrayBuffer, Set}

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

class MerSet() {

  val core_mers = Set[String]()

  val mer_counts = Map[String, Int]()

  val top_mers = Set[String]()

  var set_count = 0

  // ***************************************************************************
  def tallyMers( mers: Set[String]) {
    set_count += 1

    // Tally each of the mers.
    mers foreach {case(word) =>
      if ( mer_counts contains word )
        mer_counts( word ) += 1
      else
        mer_counts += (word -> 1)
    }  // foreach
  }  // tallyMers
 
  // ***************************************************************************
  def bestMers() {
    val half = (set_count + 1) / 2
    val qtr3 = (set_count * 3 + 1) / 4
    mer_counts foreach {case(word, count) =>
      if ( count >= half ) top_mers += word
      if ( count >= qtr3 ) core_mers += word
    }  // foreach

    println( "bestMers: (" + half + ") " + top_mers.mkString( "|" ) )
    println( "bestMers: (" + qtr3 + ") " + core_mers.mkString( "|" ) )
  }  // bestMers

  // ***************************************************************************
}  // class MerSet
