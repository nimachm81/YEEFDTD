
#ifndef FDTD_PARAMFILETRANSLATOR_H_
#define FDTD_PARAMFILETRANSLATOR_H_


#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "NumberTypes.h"
#include "ParameterExtractor.h"
#include "SingleGridParameterExtractor.h"
#include "GridCollectionParameterExtractor.h"
#include "YeeGrid.h"
#include "YeeGridCollection.h"


class ParamFileTranslator {
    public:
    ParamFileTranslator(const std::string filename);

    void Translate();
    void TranslateSingleGrid(boost::property_tree::ptree node);
    void TranslateGridCollection(boost::property_tree::ptree node);

    void SetSingleGridDimensions(YeeGrid3D& yee, SingleGridParameterExtractor& singleGridRoot);
    void SetSingleGridGridArrays(YeeGrid3D& yee, SingleGridParameterExtractor& singleGridRoot);
    void SetSingleGridGeometries(YeeGrid3D& yee, SingleGridParameterExtractor& singleGridRoot);
    void SetSingleGridParticleEmitters(YeeGrid3D& yee, SingleGridParameterExtractor& singleGridRoot);
    void SetSingleGridGirdArrayManipulatorUpdaters(YeeGrid3D& yee, SingleGridParameterExtractor& singleGridRoot);
    void SetSingleGridGridArrayManipulators(YeeGrid3D& yee, SingleGridParameterExtractor& singleGridRoot);
    void SetSingleGridUpddateInstructions(YeeGrid3D& yee,
                                          SingleGridParameterExtractor& singleGridRoot,
                                          std::map<std::string, YeeGrid3D*>& gridsMap       // maps grid names to grids
                                                                                           // in a grid-collection simulation
                                                                                           // in a single grid simulation it
                                                                                           // is not used
                                          );
    void SetSingleGridUpdateSequences(YeeGrid3D& yee, SingleGridParameterExtractor& singleGridRoot);
    void SetSingleGridGridViews(YeeGrid3D& yee, SingleGridParameterExtractor& singleGridRoot);
    void SetAndRunSingleGridRunSequencs(YeeGrid3D& yee, SingleGridParameterExtractor& singleGridRoot);
    void SetAndRunGridCollectionRunSequencs(YeeGridCollection& gridCollection,
                                            std::map<std::string, YeeGrid3D*>& gridsMap,        // grids names --> grids data
                                            std::map<std::string, std::size_t>& gridsInds,      // grids names --> grids indices
                                            GridCollectionParameterExtractor& gridCollectionRoot);

    private:
    std::string filename;

};

#endif // FDTD_PARAMFILETRANSLATOR_H_
