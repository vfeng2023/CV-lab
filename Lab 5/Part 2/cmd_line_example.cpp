#include <iostream>
#include <string>
#include <vector>

int main(int argc, char* argv[])
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << "--destination DESTINATION SOURCE" << std::endl;
        return 1;
    }
    std::vector <std::string> sources;
    std::string destination;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--destination") {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                destination = argv[i++]; // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                  std::cerr << "--destination option requires one argument." << std::endl;
                return 1;
            }  
        } else {
            sources.push_back(argv[i]);
        }
    }
    return move(sources, destination);
}