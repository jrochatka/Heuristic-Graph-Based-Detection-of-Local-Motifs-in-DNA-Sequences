
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>

using namespace std;
int k;
int score;
void loadnucleo(const string);
void start_parameters();
void seq_filter();
void divide_into_kmers();

class Nucleotide {
public:
	char base;
	int position;
	int quality;

	Nucleotide(char b, int p, int q) {
		
		base = b;
		position = p;
		quality = q;
	}
	void print_nucleo() const{
		cout << base << " ";
		cout << position << " ";
		cout << quality << endl;
	}
};

class Kmer {
public:
	vector <char> oligo;
	vector <int> positions;
	int seq_id;

	Kmer(vector <char>& o, vector <int>& p, int q) {
		oligo = o;
		positions = p;
		seq_id = q;
	}
};

vector <Nucleotide*> sequence_I;
vector <Nucleotide*> sequence_II;
vector <Nucleotide*> sequence_III;
vector <Nucleotide*> sequence_IV;
vector <Nucleotide*> sequence_V;

vector <Nucleotide*> clean_seq_I;
vector <Nucleotide*> clean_seq_II;
vector <Nucleotide*> clean_seq_III;
vector <Nucleotide*> clean_seq_IV;
vector <Nucleotide*> clean_seq_V;


vector<vector <Kmer*>> kmer_matrix;
vector <Kmer*> spectrum;

void load_nucleo(const string& filename) {
	fstream fastaFile(string("data/") + filename + ".fasta");
	if (!fastaFile.is_open()) {
		cout << "Nie mozna otworzyc pliku " << filename << endl;
		return;
	}

	fstream qualFile(string("data/") + filename + ".qual");
	if (!qualFile.is_open()) {
		cout << "Nie mozna otworzyc pliku " << filename << endl;
		return;
	}

	int vector_id = 0; 
	int position = 1;

	string fline, qline;
	while (getline(fastaFile, fline) && getline(qualFile, qline)) {
	
		istringstream f_iss(fline);
		istringstream q_iss(qline);

		char nucleo;
		int quality;
		
		if (fline[0] == '>') {
			vector_id++;
			position = 1;
			continue;
		}
			while (f_iss >> nucleo && q_iss >> quality) {

					Nucleotide* nuc = new Nucleotide(nucleo, position, quality);

					switch (vector_id)
					{
					case 1:
						sequence_I.push_back(nuc);
						position++;

						break;
					case 2:
						sequence_II.push_back(nuc);
						position++;

						break;
					case 3: 
						sequence_III.push_back(nuc);
						position++;

						break;
					case 4:
						sequence_IV.push_back(nuc);
						position++;
						break;
					case 5:
						sequence_V.push_back(nuc);
						position++;
						break;
					}
			}
	}

	fastaFile.close();
	qualFile.close();
}

void start_parameters() {

	cout << "Podaj dlugosc podciagow (k):" << endl;
	cin >> k;
	cout << endl;
	cout << "Podaj najnizsza dopuszczalna wiarygodnosc (0-40):" << endl;
	cin >> score;
	cout << "Podaj nazwe pliku (instancji, bez rozszerzenia):" << endl;
	string filename;
	cin >> filename;

	load_nucleo(filename);
}



void seq_filter() {
	for (int i = 0; i < sequence_I.size(); i++) {
		if (sequence_I[i]->quality >= score) {
			clean_seq_I.push_back(sequence_I[i]);
		}
	}
	for (int i = 0; i < sequence_II.size(); i++) {
		if (sequence_II[i]->quality >= score) {
			clean_seq_II.push_back(sequence_II[i]);
		}
	}
	for (int i = 0; i < sequence_III.size(); i++) {
		if (sequence_III[i]->quality >= score) {
			clean_seq_III.push_back(sequence_III[i]);
		}
	}
	for (int i = 0; i < sequence_IV.size(); i++) {
		if (sequence_IV[i]->quality >= score) {
			clean_seq_IV.push_back(sequence_IV[i]);
		}
	}
	for (int i = 0; i < sequence_V.size(); i++) {
		if (sequence_V[i]->quality >= score) {
			clean_seq_V.push_back(sequence_V[i]);
		}
	}
}

void divide_into_kmers(){

	if (!(clean_seq_I.size() < k)) {
		for (int i = 0; i <= (clean_seq_I.size() - k); i++) {

			vector <char>* tmp_seq = new vector <char>();
			vector <int>* tmp_position = new vector <int>();
			int seq_id = 1;

			for (int j = 0; j < k; j++) {

				tmp_seq->push_back(clean_seq_I[i + j]->base);
				tmp_position->push_back(clean_seq_I[i + j]->position);
			}

			Kmer* k = new Kmer(*tmp_seq, *tmp_position, seq_id);

			spectrum.push_back(k);

			delete tmp_seq;
			delete tmp_position;
		}
	}

	if (!(clean_seq_II.size() < k)) {
		for (int i = 0; i <= (clean_seq_II.size() - k); i++) {

			vector <char>* tmp_seq = new vector <char>();
			vector <int>* tmp_position = new vector <int>();
			int seq_id = 2;

			for (int j = 0; j < k; j++) {

				tmp_seq->push_back(clean_seq_II[i + j]->base);
				tmp_position->push_back(clean_seq_II[i + j]->position);
			}

			Kmer* k = new Kmer(*tmp_seq, *tmp_position, seq_id);

			spectrum.push_back(k);

			delete tmp_seq;
			delete tmp_position;
		}
	}

	if (!(clean_seq_III.size() < k)) {
		for (int i = 0; i <= (clean_seq_III.size() - k); i++) {

			vector <char>* tmp_seq = new vector <char>();
			vector <int>* tmp_position = new vector <int>();
			int seq_id = 3;

			for (int j = 0; j < k; j++) {

				tmp_seq->push_back(clean_seq_III[i + j]->base);
				tmp_position->push_back(clean_seq_III[i + j]->position);
			}

			Kmer* k = new Kmer(*tmp_seq, *tmp_position, seq_id);

			spectrum.push_back(k);

			delete tmp_seq;
			delete tmp_position;
		}
	}

	if (!(clean_seq_IV.size() < k)) {
		for (int i = 0; i <= (clean_seq_IV.size() - k); i++) {

			vector <char>* tmp_seq = new vector <char>();
			vector <int>* tmp_position = new vector <int>();
			int seq_id = 4;

			for (int j = 0; j < k; j++) {

				tmp_seq->push_back(clean_seq_IV[i + j]->base);
				tmp_position->push_back(clean_seq_IV[i + j]->position);
			}

			Kmer* k = new Kmer(*tmp_seq, *tmp_position, seq_id);

			spectrum.push_back(k);

			delete tmp_seq;
			delete tmp_position;
		}
	}

	if (!(clean_seq_V.size() < k)) {
		for (int i = 0; i <= (clean_seq_V.size() - k); i++) {

			vector <char>* tmp_seq = new vector <char>();
			vector <int>* tmp_position = new vector <int>();
			int seq_id = 5;

			for (int j = 0; j < k; j++) {

				tmp_seq->push_back(clean_seq_V[i + j]->base);
				tmp_position->push_back(clean_seq_V[i + j]->position);
			}

			Kmer* k = new Kmer(*tmp_seq, *tmp_position, seq_id);

			spectrum.push_back(k);

			delete tmp_seq;
			delete tmp_position;
		}
	}
}

vector<vector <int>> fill_in_matrix() {
	vector<vector <int>> graph_matrix(spectrum.size(), vector <int>(spectrum.size(), 0));
	for (int i = 0; i < spectrum.size(); i++) {
		for (int j = 0; j < spectrum.size(); j++) {

			if ((spectrum[i]->oligo == spectrum[j]->oligo) && (spectrum[i]->seq_id != spectrum[j]->seq_id && (((abs(spectrum[i]->positions[0] - spectrum[j]->positions[0])) < (k * 10))))) {
				graph_matrix[i][j] = 1;
			}
		}
	}

	return graph_matrix;
}

vector <vector<Kmer*>> possible_solutions_list;

void possible_star_structure(vector<vector <int>> graph_matrix) {
	int found = 0;


	for (int i = 0; i < spectrum.size(); i++) {
		int sum = 0;

		for (int j = 0; j < spectrum.size(); j++) {
			sum = sum + graph_matrix[i][j];
			if (sum >= 4) {
				found++;

				vector <Kmer*> tmp_solution;

				tmp_solution.push_back(spectrum[i]);

				for (int d = 0; d < spectrum.size(); d++) {
					if (graph_matrix[i][d] == 1) {
						tmp_solution.push_back(spectrum[d]);
					}
				}
				
				possible_solutions_list.push_back(tmp_solution);

			}
	
		}

	}

	cout << "Znaleziono tyle mozliwych rozwiazan: " << found << endl;
	}

vector <Kmer*> final_solution;
void is_from_different_seq() {
	

	for (int j = 0; j < possible_solutions_list.size(); j++) {

		bool seqI = false;
		bool seqII = false;
		bool seqIII = false;
		bool seqIV = false;
		bool seqV = false;

		bool seqI_final = false;
		bool seqII_final = false;
		bool seqIII_final = false;
		bool seqIV_final = false;
		bool seqV_final = false;


		for (int h = 0; h < possible_solutions_list[j].size(); h++) {

			int seq_id = possible_solutions_list[j][h]->seq_id;

			if (seq_id == 1) seqI = true;
			if (seq_id == 2) seqII = true;
			if (seq_id == 3) seqIII = true;
			if (seq_id == 4) seqIV = true;
			if (seq_id == 5) seqV = true;

		}

		if (seqI && seqII && seqIII && seqIV && seqV) {

			for (int v = 0; v < possible_solutions_list[j].size(); v++) {

				if ((possible_solutions_list[j][v]->seq_id == 1) && (seqI_final == false)) {

					final_solution.push_back(possible_solutions_list[j][v]);
					seqI_final = true;
					continue;
				}

				if ((possible_solutions_list[j][v]->seq_id == 2) && (seqII_final == false)) {

					final_solution.push_back(possible_solutions_list[j][v]);
					seqII_final = true;
					continue;
				}

				if ((possible_solutions_list[j][v]->seq_id == 3) && (seqIII_final == false)) {

					final_solution.push_back(possible_solutions_list[j][v]);
					seqIII_final = true;
					continue;
				}

				if ((possible_solutions_list[j][v]->seq_id == 4) && (seqIV_final == false)) {

					final_solution.push_back(possible_solutions_list[j][v]);
					seqIV_final = true;
					continue;
				}

				if ((possible_solutions_list[j][v]->seq_id == 5) && (seqV_final == false)) {

					final_solution.push_back(possible_solutions_list[j][v]);
					seqV_final = true;
					continue;
				}
			}
			return;
		}
		
	}
}

void static result_display() {

	if (!final_solution.empty()) {
		cout << "Rozwiazanie: " << endl;
		for (int i = 0; i < final_solution.size(); i++) {

			cout << " Numer sekwencji: " << final_solution[i]->seq_id << " Pozycja k-meru w sekwencji: " << final_solution[i]->positions[0] << " Sekwencja nukleotydowa k-meru: ";

			for (int j = 0; j < final_solution[i]->oligo.size(); j++) {
				cout << final_solution[i]->oligo[j];
			}
			cout << endl;
		}
	}
	else { cout << "Nie znaleziono rozwiązania :(" << endl; }
}


int main()
{
	start_parameters();
	seq_filter();
	divide_into_kmers();
	fill_in_matrix();
	possible_star_structure(fill_in_matrix());
	is_from_different_seq();
	result_display();
}

