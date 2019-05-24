
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
class AminoBase( val name: String, val chain_id: Char, val res_seq: String, val icode: Char ) {

  val translate: Map[String, Char] = Map(
    "ALA" -> 'A', "ARG" -> 'R', "ASN" -> 'N', "ASP" -> 'D', "ASX" -> 'B', "CYS" -> 'C',
    "GLU" -> 'E', "GLN" -> 'Q', "GLX" -> 'Z', "GLY" -> 'G', "HIS" -> 'H', "ILE" -> 'I',
    "LEU" -> 'L', "LYS" -> 'K', "MET" -> 'M', "PHE" -> 'F', "PRO" -> 'P', "SER" -> 'S',
    "THR" -> 'T', "TRP" -> 'W', "TYR" -> 'Y', "VAL" -> 'V')

  val letter = if ( translate contains name ) translate( name ) else 'X'

  // **************************************************************************
  override def toString: String = {
    name + "|" + letter + "|" + chain_id + "|" + res_seq + "." + icode + "|"
  }  // toString

  // **************************************************************************
}  // class AminoBase
