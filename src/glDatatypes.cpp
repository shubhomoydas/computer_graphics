#include "glDatatypes.h"
#include <iostream>
#include <fstream>
#include <string>

namespace smd {

	/**
	 * Vector
	 */
	bool Vector::parse(std::string str) {
		if (!str.empty()) {
			std::string nstr = trim(str);
			std::vector<std::string> tokens;
			tokenize(nstr, tokens, ",");
			if (tokens.size() >= 3) {
				x = std::atof(tokens[0].c_str());
				y = std::atof(tokens[1].c_str());
				z = std::atof(tokens[2].c_str());
				w = 0;
				if (tokens.size() > 3) w = std::atof(tokens[3].c_str());
				return true;
			} else {
				//std::cout << "Invalid vector: " << str << std::endl;
			}
		}
		return false;
	}

	/**
	 * Color
	 */
	bool Color::parse(std::string str) {
		if (!str.empty()) {
			std::string nstr = trim(str);
			std::vector<std::string> tokens;
			tokenize(nstr, tokens, ",");
			if (tokens.size() >= 3) {
				r = std::atof(tokens[0].c_str());
				g = std::atof(tokens[1].c_str());
				b = std::atof(tokens[2].c_str());
				s = 0;
				if (tokens.size() > 3) s = std::atof(tokens[3].c_str());
				return true;
			} else {
				//std::cout << "Invalid color: " << str << std::endl;
			}
		}
		return false;
	}

	/**
	 *  MatrixStack
	 */
	void MatrixStack::initialize() {
		clear();
		top = new Matrix(4);
		top->setIdentity();
		mats.push_back(top);
	}
	void MatrixStack::push() {
		mats.push_back(top->clone());
		top = mats.back();
	}
	Matrix *MatrixStack::pop() {
		if (mats.size() < 2) {
			std::cout << "MatrixStack: Error! Pop attempted with only " << 
				mats.size() << " element(s) on stack." << std::endl;
			exit(-1);
		}

		// delete current top
		Matrix *m = mats.back();
		mats.pop_back();
		delete m;

		top = mats.back();
		return top;
	}
	int MatrixStack::size() {
		return mats.size();
	}
	void MatrixStack::clear() {
		top = 0;
		while (!mats.empty()) {
			Matrix *m = mats.front();
			mats.pop_front();
			delete m;
		}
	}

	MatrixStack::~MatrixStack() {
		clear();
	}

	std::ostream& operator<< (std::ostream& out, const Vector& vector) {
		vector.print(out);
		return out;
	}

	std::ostream& operator<< (std::ostream& out, const Color& color) {
		return color.print(out);
	}

	std::ostream& operator<< (std::ostream& out, const Vertex& vertex) {
		return vertex.print(out);
	}

	std::string trim(const std::string& str, const std::string& chars)
	{
		std::string::size_type pos1 = str.find_first_not_of(' ');
		std::string::size_type pos2 = str.find_last_not_of(' ');
		std::string newstr = str.substr(pos1 == std::string::npos ? 0 : pos1,
		pos2 == std::string::npos ? str.length() - 1 : pos2 - pos1 + 1);
		return newstr;
	}

	void tokenize(const std::string& str,
						std::vector<std::string>& tokens,
						const std::string& delimiters) {
		// Skip delimiters at beginning.
		std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
		// Find first "non-delimiter".
		std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

		while (std::string::npos != pos || std::string::npos != lastPos)
		{
			// Found a token, add it to the vector.
			std::string token = str.substr(lastPos, pos - lastPos);
			trim(token);
			tokens.push_back(token);
			// Skip delimiters.  Note the "not_of"
			lastPos = str.find_first_not_of(delimiters, pos);
			// Find next "non-delimiter"
			pos = str.find_first_of(delimiters, lastPos);
		}

	}

	void saveData(const Matrix& mat, const std::string& filepath) {
		if (filepath.empty()) {
			std::cerr << "Output filename missing..." << std::endl;
			return;
		}
		double **a = mat.data;
		std::ofstream outfile(filepath.c_str());
		for (int i = 0; i < mat.rows; i++) {
			for (int j = 0; j < mat.cols; j++) {
				if (j > 0) outfile << ",";
				outfile << a[i][j];
			}
			outfile << std::endl;
		}
		outfile.close();
	}

	Matrix* populateMatFromStringVectorOfFloats(const std::vector< std::vector<std::string>* >& recs, const int skipleftcols, const int skiprightcols) {
		int cols = 0;
		if (recs.size() > 0) {
			cols = recs[0]->size();
		}
		Matrix* mat = NULL;

		mat = new Matrix(recs.size(), cols-skipleftcols-skiprightcols);
		double **a = mat->data;
		int i = 0;
		for ( std::vector< std::vector<std::string>* >::const_iterator it=recs.begin() ; it < recs.end(); it++, i++ ) {
			std::vector<std::string>* tokens = *it;
			int j = 0;
			for ( std::vector<std::string>::iterator itt=tokens->begin() ; itt < tokens->end(); itt++, j++ ) {
				if (j >= cols-skiprightcols) {
					break;
				}
				if (j >= skipleftcols) {
					a[i][j-skipleftcols] = atof((*itt).c_str());
				}
			}
		}
		return mat;
	}

	Matrix* loadCSV(const std::string& filename, const int skiptoplines, const int skipleftcols, const int skiprightcols) {
		Matrix* mat = NULL;
		std::ifstream ifile (filename.c_str());
		int skip = skiptoplines;
		if (ifile.is_open()) {
			std::vector<std::vector<std::string>* >* recs = new std::vector< std::vector<std::string>* >;
			int reccount = 0;
			//cout << "Opened the file..." << endl;
			while ( ifile.good() ) {
				std::string line;
				getline (ifile,line);
				if (skip > 0) {
					skip--;
					continue;
				}
				trim(line);
				if (line.length() == 0) continue;
				reccount++;
				std::vector<std::string>* tokens = new std::vector<std::string>;
				tokenize(line, *tokens, ",");
				recs->push_back(tokens);
				if (reccount == 1) {
					// get the number of columns
					//cout << "Columns: " << tokens->size() << "; line length: " << line.length() << endl;
				}
			}
			ifile.close();
			//cout << "Total rows read: " << reccount << endl;
			mat = populateMatFromStringVectorOfFloats(*recs, skipleftcols, skiprightcols);
			release_vector_elements< std::vector<std::string>* >(recs);
			delete recs;
		} else std::cout << "Unable to open file " << filename << std::endl;
		return mat;
	}

} // namespace smd
