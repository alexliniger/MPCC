#!/bin/bash
echo "--- Preparing Docker environment ---"
mkdir -p .gitman_cache

# Define the command to run inside the container
# We use a marker file check to be as fast as possible
CONTAINER_CMD="
cd C++
if [ ! -f 'External/json/CMakeLists.txt' ]; then
    echo '--- Initializing dependencies with gitman (this may take a moment) ---'
    sudo -u devuser gitman install --force
else
    echo '--- Dependencies detected, skipping initialization ---'
fi
cd ..

if [ ! -f '.git/hooks/pre-commit' ]; then
    echo '--- Installing pre-commit hooks ---'
    sudo -u devuser pre-commit install --config C++/.pre-commit-config.yaml > /dev/null 2>&1
fi

echo '--- Interactive shell ready ---'
echo 'Hint: You are in the C++ subdirectory. Repository root is at ..'
echo ''
cd C++
exec sudo -u devuser /bin/bash
"

# Handle SSH_AUTH_SOCK for Mac/Linux
SSH_AUTH_ARGS=""
if [ -n "$SSH_AUTH_SOCK" ]; then
    SSH_AUTH_ARGS="-v $SSH_AUTH_SOCK:/ssh-agent -e SSH_AUTH_SOCK=/ssh-agent"
fi

echo "--- Launching container ---"
docker run -it --rm \
    -p 8080:8080 \
    $SSH_AUTH_ARGS \
    -v "$(pwd)/..:/app/MPCC" \
    -v "$(pwd)/.gitman_cache:/home/devuser/.gitcache" \
    mpcc /bin/bash -c "$CONTAINER_CMD"
