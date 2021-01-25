#pragma once

enum LayerNodeFlagMask {
    VALID = (1<<0),
    IS_GHOST = (1<<1),
    GHOST_FROM_COARSE = (1<<2), // only valid when IS_GHOST flag is on.
};