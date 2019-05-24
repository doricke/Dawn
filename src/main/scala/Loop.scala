
import scala.collection.mutable.{ArrayBuffer, Set}

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

class Loop( val start: Int, val end: Int, val gaps: ArrayBuffer[Gap], val inserts: ArrayBuffer[Insert] ) {

// ****************************************************************************

  var after: Int = 0		// select insertion position

  var size: Int = 0		// insertion size

  val names = Set[String]()	// sequences spanned by this loop

// ****************************************************************************
  def findGapSite() {
    val sites = ArrayBuffer[Int]()
    gaps foreach { case(gap) => sites += gap.start }

    val sorted = sites.sorted
    // println( "findGapSize: sites.size: " + sites.size + ", sorted: " + sorted.mkString( "|" ) )
    after = sorted( sites.size/2 )
  }  // findGapSite

// ****************************************************************************
  def findSite() {
    val sites = ArrayBuffer[Int]()
    inserts foreach { case(ins) => sites += ins.after }
    // gaps foreach { case(gap) => sites += (gap.start-1) }

    val sorted = sites.sorted
    // println( "findSize: sites.size: " + sites.size + ", sorted: " + sorted.mkString( "|" ) )
    after = sorted( sites.size/2 )
    findSize()
  }  // findSite

// ****************************************************************************
  def findSize() {
    inserts foreach { case(ins) => if (ins.seq.size > size) size = ins.seq.size }
  }  // findSize

// ****************************************************************************

}  // class Loop
