
#ifndef FDTD_PARAMFILETRANSLATOR_H_
#define FDTD_PARAMFILETRANSLATOR_H_


#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "NumberTypes.h"
#include "ParameterExtractor.h"
#include "SingleGridParameterExtractor.h"


class ParamFileTranslator {
    public:
    ParamFileTranslator(const std::string filename);

    void Translate();
    void TranslateSingleGrid(boost::property_tree::ptree node);

    private:
    std::string filename;

};

#endif // FDTD_PARAMFILETRANSLATOR_H_
