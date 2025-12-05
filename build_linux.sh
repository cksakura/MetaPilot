#!/bin/bash

###############################################################################
# MetaPilot Linux Build Script
# 
# This script builds the MetaPilot web application for Linux deployment.
# It creates a production-ready JAR file and optionally packages it with
# all dependencies.
#
# Usage:
#   ./build_linux.sh              # Build web application
#   ./build_linux.sh clean        # Clean before building
#   ./build_linux.sh all          # Build all modules
#
###############################################################################

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
PROJECT_NAME="MetaPilot"
VERSION="1.0.0-beta.3"
BUILD_PROFILE="linux-web"
MAVEN_OPTS="-Xmx2g"

# Functions
print_header() {
    echo -e "${BLUE}================================================${NC}"
    echo -e "${BLUE}  $PROJECT_NAME Linux Build${NC}"
    echo -e "${BLUE}  Version: $VERSION${NC}"
    echo -e "${BLUE}================================================${NC}"
    echo ""
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ $1${NC}"
}

print_info() {
    echo -e "${YELLOW}→ $1${NC}"
}

check_requirements() {
    print_info "Checking requirements..."
    
    # Check Java version
    if ! command -v java &> /dev/null; then
        print_error "Java is not installed or not in PATH"
        exit 1
    fi
    
    JAVA_VERSION=$(java -version 2>&1 | grep -i version | head -1)
    print_success "Found Java: $JAVA_VERSION"
    
    # Check Maven
    if ! command -v mvn &> /dev/null; then
        print_error "Maven is not installed or not in PATH"
        exit 1
    fi
    
    MAVEN_VERSION=$(mvn -version | head -1)
    print_success "Found Maven: $MAVEN_VERSION"
    
    echo ""
}

clean_build() {
    print_info "Cleaning previous builds..."
    mvn clean -q
    print_success "Clean complete"
    echo ""
}

build_web_app() {
    print_info "Building MetaPilot Web Application..."
    
    mvn clean package -P$BUILD_PROFILE -DskipTests $MAVEN_OPTS
    
    if [ $? -eq 0 ]; then
        print_success "Build completed successfully"
    else
        print_error "Build failed"
        exit 1
    fi
    echo ""
}

build_all_modules() {
    print_info "Building all modules..."
    
    mvn clean package -Pall -DskipTests $MAVEN_OPTS
    
    if [ $? -eq 0 ]; then
        print_success "Build completed successfully"
    else
        print_error "Build failed"
        exit 1
    fi
    echo ""
}

show_build_info() {
    print_info "Build artifacts created:"
    echo ""
    
    if [ -f "metapilot-web/target/metapilot-web-$VERSION.jar" ]; then
        JAR_FILE="metapilot-web/target/metapilot-web-$VERSION.jar"
        FILE_SIZE=$(du -h "$JAR_FILE" | cut -f1)
        
        echo -e "  ${GREEN}✓${NC} MetaPilot Web JAR"
        echo -e "     Location: $JAR_FILE"
        echo -e "     Size: $FILE_SIZE"
        echo ""
        echo -e "  ${BLUE}To run:${NC}"
        echo -e "     java -jar $JAR_FILE"
        echo ""
    else
        print_error "JAR file not found"
    fi
    
    # Check Windows GUI JAR if building all
    if [ -f "metapilot-windows-gui/target/metapilot-windows-gui-$VERSION-shaded.jar" ]; then
        GUI_JAR="metapilot-windows-gui/target/metapilot-windows-gui-$VERSION-shaded.jar"
        GUI_SIZE=$(du -h "$GUI_JAR" | cut -f1)
        
        echo -e "  ${GREEN}✓${NC} MetaPilot Windows GUI JAR"
        echo -e "     Location: $GUI_JAR"
        echo -e "     Size: $GUI_SIZE"
        echo ""
    fi
}

# Main execution
main() {
    print_header
    check_requirements
    
    # Parse arguments
    case "${1:-build}" in
        clean)
            clean_build
            build_web_app
            ;;
        all)
            build_all_modules
            ;;
        build)
            build_web_app
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [clean|all|build]"
            exit 1
            ;;
    esac
    
    show_build_info
    
    echo -e "${GREEN}================================================${NC}"
    echo -e "${GREEN}  Build Complete!${NC}"
    echo -e "${GREEN}================================================${NC}"
}

# Run main function
main "$@"

