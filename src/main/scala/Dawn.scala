
//  This is the main class for Dawn class.
//
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

// ******************************************************************************
object Dawn {

  //****************************************************************************
  // MSA program.
  def main(args: Array[String]) {

    val seqs = if ( args.length >= 1 ) args(0) else "ha"

    val targets = if ( args.length >= 2 ) args(1) else ""

    val pdbs = if ( args.length >= 3 ) args(2) else ""

    val dealer: MsaSetup = new MsaSetup( seqs, targets, pdbs )
    dealer.process()
  }  // main

}  // object Dawn
