#include <stdio.h>
#include <unistd.h>
#include <string>
#include <vector>
#include <map>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


typedef uint32_t Color;
typedef uint32_t TileIndex;

struct Parameters
{
    uint16_t tileWidth, tileHeight;
    std::string filename;
};

struct Tile
{
    std::vector<Color> bitmap;
    //TileIndex tileIndex;
    
    uint64_t ComputeChecksum() const
    {
        uint64_t res = 0;
        for(const auto color : bitmap)
        {
            res += color;
        }
        return res;
    }
    
    bool operator == (const Tile& other) const
    {
        assert(bitmap.size() == other.bitmap.size());
        for(size_t i = 0; i < bitmap.size(); i++)
        {
            if (bitmap[i] != other.bitmap[i])
            {
                return false;
            }
        }
        return true;
    }
    
    void BuildFromMap(Color* colors, uint32_t offset, uint32_t pitch, const Parameters& parameters)
    {
        const auto w = parameters.tileWidth;
        const auto h = parameters.tileHeight;
        bitmap.resize(h * w);
        auto* p = bitmap.data();
        for (size_t y = 0; y < h; y++)
        {
            for (size_t x = 0; x < w; x++)
            {
                *p++ = colors[offset + y * pitch +x];
            }
        }
    }
};

std::vector<Tile> tiles;
std::map<uint64_t, std::vector<TileIndex> > uniqueTiles;

TileIndex GetTileIndex(Color* colors, uint32_t offset, uint32_t pitch, const Parameters& parameters)
{
    Tile tile;
    tile.BuildFromMap(colors, offset, pitch, parameters);
    auto checkSum = tile.ComputeChecksum();
    
    auto iter = uniqueTiles.find(checkSum);
    if (iter != uniqueTiles.end())
    {
        auto& existingTileIndices = iter->second;
        for (const auto tileIndex : existingTileIndices)
        {
            if (tiles[tileIndex] == tile)
            {
                // tile found
                return tileIndex;
            }
        }
    }
    // not found, add new
    TileIndex res = tiles.size();
    tiles.push_back(tile);
    uniqueTiles[checkSum].push_back(res);
}

bool ParseParameters(int argc, char** argv, Parameters& parameters)
{
    parameters.tileWidth = 16;
    parameters.tileHeight = 16;
    
    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd)) != NULL)
    {
        printf("Current working dir: %s\n", cwd);
    }
    
    parameters.filename = std::string(cwd) + "/../../maps/Zelda3LightOverworldBG_masked.png";
    return true;
}

int main(int argc, char** argv)
{
    printf("WFCTiles\n");
    
    Parameters parameters;
    if (!ParseParameters(argc, argv, parameters))
    {
        return 1;
    }

    int x,y,n;
    unsigned char *data = stbi_load(parameters.filename.c_str(), &x, &y, &n, 4);
    if (!data)
    {
        return 1;
    }
    stbi_image_free(data);
    return 0;
}
