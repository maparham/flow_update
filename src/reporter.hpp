#ifndef REPORTER_HPP_
#define REPORTER_HPP_

#define RESULT_DIR string("data/output/")

template<typename Out>
void split(const std::string &s, char delim, Out result) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		*(result++) = item;
	}
}

std::vector<std::string> split(const std::string &s, char delim) {
	std::vector < std::string > elems;
	split(s, delim, std::back_inserter(elems));
	return elems;
}

struct Reporter {
	string name;
	static int counter;
	int index;
	ofstream of, all;

	static string resultDir() {
		return RESULT_DIR;
	}

	// constructor
	Reporter(const string path) {
		// prepare result file
		vector < string > tokens = split(path, '/');
		string dir = resultDir();
		string resultPath = RESULT_DIR + tokens.back() + ".txt";
		name = tokens.back();
		printf("Result Path: %s\n", resultPath.c_str());
		of.open(resultPath);
		all.open(RESULT_DIR + "all.txt", std::ios_base::app);
		if (!of.is_open()) {
			printf("file not opened: %s", resultPath.c_str());
			exit(1);
		}
		index = counter++;
	}
	template<typename T>
	Reporter& operator<<(T input) {
		of << input;
		all << input;
		return *this;
	}
	~Reporter() {
		of.close();
	}
};
int Reporter::counter = 0;

#endif /* REPORTER_HPP_ */
