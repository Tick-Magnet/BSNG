// Usage ./generator random_points.bin random_points.txt 20 1 20
// Compile: g++ -o generator generator.cpp

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <string>

struct Point {
    float x;
    float y;
};

void generate_random_points(const std::string& bin_file_name, const std::string& txt_file_name, int num_points, int min_value, int max_value) {
    std::ofstream bin_file(bin_file_name, std::ios::binary);
    std::ofstream txt_file(txt_file_name);

    if (!bin_file || !txt_file) {
        std::cerr << "Error opening the files." << std::endl;
        return;
    }

    // Write the number of points (N) and number of dimensions (D) as 4-byte integers
    bin_file.write(reinterpret_cast<const char*>(&num_points), sizeof(int)); // N
    int num_dimensions = 2; // Two dimensions X/Y
    bin_file.write(reinterpret_cast<const char*>(&num_dimensions), sizeof(int)); // D

    // Random number generator setup
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(min_value, max_value);

    // Vector to store the points
    std::vector<Point> points;

    // Generate random points and write to files
    for (int i = 0; i < num_points; ++i) {
        Point point;
        point.x = static_cast<float>(dist(gen));
        point.y = static_cast<float>(dist(gen));
        points.push_back(point);

        // Write the points' coordinates to the text file
        txt_file << "Point " << i + 1 << ": " << point.x << ", " << point.y << std::endl;
    }

    // Write the points' coordinates to the binary file
    bin_file.write(reinterpret_cast<const char*>(points.data()), points.size() * sizeof(Point));

    bin_file.close();
    txt_file.close();
}

int main(int argc, char* argv[]) {
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " <bin_file_name> <txt_file_name> <num_points> <min_value> <max_value>" << std::endl;
        return 1;
    }

    // Command-line arguments
    std::string bin_file_name = argv[1];
    std::string txt_file_name = argv[2];
    int num_points = std::stoi(argv[3]);
    int min_value = std::stoi(argv[4]);
    int max_value = std::stoi(argv[5]);

    if (num_points <= 0 || min_value >= max_value) {
        std::cerr << "Invalid arguments. Please provide a valid number of points and a valid range." << std::endl;
        return 1;
    }

    // Generate and save the .bin and .txt files
    generate_random_points(bin_file_name, txt_file_name, num_points, min_value, max_value);

    return 0;
}