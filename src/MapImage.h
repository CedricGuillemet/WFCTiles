#pragma once
#include <vector>
#include "Parameters.h"

typedef uint32_t Color;

void Blit(const Color* tileInput, Color* output, uint32_t offset, uint32_t pitch, uint16_t tileWidth, uint16_t tileHeight)
{
    for (size_t y = 0; y < tileHeight; y++)
    {
        for (size_t x = 0; x < tileWidth; x++)
        {
            Color* destination = output + pitch * y + offset + x;
            *destination = *tileInput++;
        }
    }
}

void SaveMapImage(const Parameters& parameters, int mapIndex)
{
    const auto w = parameters.mapWidth;
    const auto h = parameters.mapHeight;
    std::vector<Color> image(w * parameters.tileWidth * h * parameters.tileHeight);

    for (size_t ty = 0; ty < w; ty++)
    {
        for (size_t tx = 0; tx < h; tx++)
        {
            Coord coord{ int32_t(tx), int32_t(ty) };
            const auto tileIndex = GetTileAtIndex(coord);
            const auto& tile = tiles[tileIndex];
            const uint32_t offset = uint32_t(ty * parameters.tileHeight * w * parameters.tileWidth + tx * parameters.tileWidth);
            Blit(tile.bitmap.data(), image.data(), offset, w * parameters.tileWidth, parameters.tileWidth, parameters.tileHeight);
        }
    }

    // save image
    const auto imageWidth = w * parameters.tileWidth;
    char mapFilename[512];
    sprintf(mapFilename, "MapGen_seed%d_%d_%dx%d.png", parameters.seed, mapIndex, int(parameters.mapWidth), int(parameters.mapHeight));
#ifdef _MSC_VER
    int res = stbi_write_png(mapFilename, imageWidth, h * parameters.tileHeight, 4,
        image.data(), imageWidth * 4);
#else
    int res = stbi_write_png("/Users/cedricguillemet/dev/WFCTiles/res.png", imageWidth, mHeight * parameters.tileHeight, 4,
        image.data(), imageWidth * 4);
#endif
}