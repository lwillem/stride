#############################################################################
# COMMAND LINE INTERFACE TOOLS
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright (C) 2019 lwillem, SIMID, UNIVERSITY OF ANTWERP, BELGIUM
#############################################################################

#' @title Convert a list into XML format (RECURSIVE FUNCTION)
#'
#' @description The counterpart of the "XML::xmlToList" function
#'
#' @param node     the XML node
#' @param sublist  the list to be stored
#'
#' @export
#' @import XML
smd_listToXML <- function(node, sublist){

  for(i in 1:length(sublist)){
    child <- newXMLNode(names(sublist)[i], parent=node);

    if (length(sublist[[i]]) > 1){
      tmp <- unlist(list(sublist[[i]]))
      smd_listToXML(child, tmp)
    }
    else{
      xmlValue(child) <- paste(sublist[[i]],collapse= ';')
    }
  }
}

#' @title Save a list as XML file with given root node
#'
#' @description The counterpart of the "XML::xmlToList" function
#'
#' @param data_list          the list to be stored
#' @param root_name          the name of the XML main root
#' @param file_name_prefix   the prefix for the file name
#' @param xml_prefix         the xml prefix (optional)
#'
#' @keywords external
#' @export
#' @import XML
smd_save_as_xml <- function(data_list,
                            root_name,
                            file_name_prefix = 'data',
                            xml_prefix = NULL){

  # setup XML doc (to add prefix)
  xml_doc <- newXMLDoc()

  # setup XML root
  root    <- newXMLNode(root_name, doc = xml_doc)

  # add list info
  smd_listToXML(root, data_list)

  # create filename
  file_name <- smd_file_path(paste0(file_name_prefix,'.xml'))

  # xml prefix
  if(!is.null(xml_prefix)){
    xml_prefix <- newXMLCommentNode(xml_prefix)
  }

  # save as XML,
  # note: if we use an XMLdoc to include prefix, the line break dissapears...
  # fix: http://r.789695.n4.nabble.com/saveXML-prefix-argument-td4678407.html
  cat(saveXML(xml_doc, indent = TRUE, prefix = xml_prefix), file = file_name)

  # return the file name
  return(file_name)
}

