
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

class Conserved() {

  // I:L:V, R:K, F:Y
  val allowed = Map( 'A' -> Set('S'), 'D' -> Set('E', 'N'), 'E' -> Set('D', 'Q'), 
      'F' -> Set('Y'), 'I' -> Set('L', 'V'), 'K' -> Set('R'), 'L' -> Set('I', 'V'), 
      'N' -> Set('D', 'Q'), 'Q' -> Set('E', 'N'), 'R' -> Set('K'), 
      'S' -> Set('A', 'T'), 'T' -> Set('S'), 'Y' -> Set('F'), 'V' -> Set('I', 'L') )

  val rules = Array(Set('I', 'L', 'V'), Set('R', 'K'), Set('S', 'T'), Set('S', 'A'), Set('E', 'D'), Set('D', 'N'), Set('Q', 'E'), Set('Q', 'N'), Set('F', 'Y'))

  // ***************************************************************************
  def conservative( observed: Map[Char, Int] ): Char = {
    if ( observed.keys.size == 0 ) 
      '.'
    else {
      var cons = '.'
      rules foreach {case(rule) =>
        val diff_set = observed.keySet -- rule
        if ( diff_set.size == 0 ) {
          cons = observed.toList.sortBy{_._2}.last._1.toLower
        }  // if
      }  // foreach
      cons
    }  // if
  }  // conservative

  // ***************************************************************************
  // This method counts up the observed nonconserved residues at this position.
  def nonconserved( aa: Char, observed: Map[Char, Int] ): Int = {
    var nc = 0
    val c_keys = allowed.keySet
    observed foreach {case(residue, count) =>
      if ( ( residue != aa ) && ( residue != '-' ) ) {
        if ( ( c_keys contains residue ) == false )
          nc += count
        else {
          if ( ( ( allowed contains aa ) == false ) || ( allowed( aa ) contains residue ) == false )
            nc += count
        }  // if
      }  // if
    }  // foreach
    nc
  }  // nonconserved

  // ***************************************************************************

}  // class Conserved

