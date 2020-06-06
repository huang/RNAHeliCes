#ifndef FILTER_ADD_HH
#define FILTER_ADD_HH

//#include "rtlib/rope.hh"

//#include "string.hh"
//#include "sequence.hh"

//#include <cassert>
#include <iostream>
#include <vector>


template<typename alphabet, typename pos_type, typename T>
inline bool exact_pairing_center_at(const Basic_Sequence<alphabet, pos_type> &seq,
    T i, T j, std::vector<std::string>& tokens)
{
    if (tokens.empty())
    {
        return (i+3 < j) && basepairing(seq, i, j) && basepairing(seq, i+1, j-1);
    }
    else
    {
	std::vector<std::string>::iterator it;
	float helix_center_float = (i+j+1)/2.0f;

	bool found = false;
	float fuzzy = 0.0f;
	for ( it=tokens.begin() ; it < tokens.end(); it++ ) {
	    fuzzy = atof((*it).c_str());
	    if ((helix_center_float >= (fuzzy-0.1f)) && (helix_center_float <= (fuzzy+0.1f))) {
		found = true;
		break;
	    }
	}
	
	return (i+3 < j) && basepairing(seq, i, j) && basepairing(seq, i+1, j-1) && found;  //((i+j+1)/2 == l);
    }
}

template<typename alphabet, typename pos_type, typename T>
inline bool pairing_center_at(const Basic_Sequence<alphabet, pos_type> &seq,
    T i, T j, std::vector<std::string>& tokens, float theta)      // modify line 1: match_str ==> str
{
    //std::cout << "line_2" << std::endl;
    if (tokens.empty())
    {
        return (i+3 < j) && basepairing(seq, i, j) && basepairing(seq, i+1, j-1);
    }
    else
    {
	std::vector<std::string>::iterator it;
	/*
		    Rope helix_index;
		    int helix_center_int;
		    helix_center_int = (i+j+1)/2;
		    if ( helix_center_int*2 > i+j+1 ) helix_center_int = helix_center_int - 1;  
		    append(helix_index, helix_center_int);
		    if ( helix_center_int*2 != i+j+1 ) append(helix_index, ".5", 2);
	  std::ostringstream outs;    // Declare an output string stream.
	  outs << helix_index;            // Convert value into a string.
	  std::string helix_index_str = outs.str();     // Get the created string from the output stream.
	  */
	float helix_center_float = (i+j+1)/2.0f;
//	std::cout << helix_center_float << std::endl;
	
// 	int helix_center_int = (i+j+1)/2;
// 	if ( helix_center_int*2 == (i+j+1) ) helix_center_int = helix_center_int - 1;  
// 	std::stringstream ss;
// 	ss << helix_center_int;
// 	if ( helix_center_int*2 != (i+j+1) ) ss << ".5";
	

	bool found = false;
	float fuzzy = 0.0f;
	for ( it=tokens.begin() ; it < tokens.end(); it++ ) {
	    fuzzy = atof((*it).c_str());
	    if ((helix_center_float >= (fuzzy-theta)) && (helix_center_float <= (fuzzy+theta))) {
// 	    if ( (ss.str() == std::to_string(value-1)) || (ss.str() == std::to_string(value-0.5)) || (ss.str() == std::to_string(value)) || (ss.str() == std::to_string(value+0.5)) || (ss.str() == std::to_string(value+1)) ) {
		found = true;
		break;
	    }
	}
	
	//std::cout << seq << std::endl; 
	//std::cout << '<' << seq.i << ", " << seq.j << '>' << std::endl;
	//std::cout << '<' << i << ", " << j << '>' << std::endl;
	//for (typename Basic_Sequence<alphabet, pos_type>::iterator i = seq.begin();i != seq.end(); ++i)
        //    std::cout << i << std::endl;  // = std::toupper(*i);
	
	return (i+3 < j) && basepairing(seq, i, j) && basepairing(seq, i+1, j-1) && found;  //((i+j+1)/2 == l);
    }
}

#endif
