#include <exception>

struct InvalidActivation : public std::exception{
    const char * what () const throw () {
        return "No Activation Information";
    }
};