#pragma once

struct Parameters
{
    uint16_t tileWidth, tileHeight;
    std::string filename;
    uint16_t mapWidth, mapHeight;
    int seed;
};


bool ParseParameters(int argc, char** argv, Parameters& parameters)
{
    parameters.tileWidth = 16;
    parameters.tileHeight = 16;

    parameters.mapWidth = 40;
    parameters.mapHeight = 10;

#ifdef _MSC_VER
    //parameters.filename = "../maps/Zelda3LightOverworldBG_masked.png";
    parameters.filename = "../maps/Zelda3LightOverworldBG.png";
#else
    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) != NULL)
    {
        printf("Current working dir: %s\n", cwd);
    }
    parameters.filename = std::string(cwd) + "/../../maps/Zelda3LightOverworldBG_masked.png";
#endif

    parameters.seed = 1337;
    return true;
}
