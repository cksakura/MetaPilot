#!/bin/bash

###############################################################################
# MetaPilot Linux Deployment Script
# 
# This script deploys the MetaPilot web application to a Linux VM.
# It handles installation, configuration, and service setup.
#
# Usage:
#   ./deploy_linux.sh                    # Deploy to default location
#   ./deploy_linux.sh /opt/metapilot    # Deploy to custom location
#
###############################################################################

set -e

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m'

# Configuration
DEPLOY_DIR="${1:-/opt/metapilot}"
JAR_FILE="metapilot-web/target/metapilot-web-1.0.0-beta.2.jar"
SERVICE_USER="metapilot"
SERVICE_NAME="metapilot-web"
CONFIG_DIR="$DEPLOY_DIR/config"
LOG_DIR="$DEPLOY_DIR/logs"

print_header() {
    echo -e "${BLUE}================================================${NC}"
    echo -e "${BLUE}  MetaPilot Linux Deployment${NC}"
    echo -e "${BLUE}================================================${NC}"
    echo ""
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_info() {
    echo -e "${YELLOW}→ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ $1${NC}"
}

check_jar_exists() {
    if [ ! -f "$JAR_FILE" ]; then
        print_error "JAR file not found: $JAR_FILE"
        echo "Please run build_linux.sh first"
        exit 1
    fi
    print_success "Found JAR file"
}

create_directories() {
    print_info "Creating deployment directories..."
    
    sudo mkdir -p "$DEPLOY_DIR"
    sudo mkdir -p "$CONFIG_DIR"
    sudo mkdir -p "$LOG_DIR"
    sudo mkdir -p "$DEPLOY_DIR/data"
    
    print_success "Directories created"
}

copy_files() {
    print_info "Copying files..."
    
    sudo cp "$JAR_FILE" "$DEPLOY_DIR/"
    sudo cp -r config/* "$CONFIG_DIR/" 2>/dev/null || true
    
    print_success "Files copied"
}

set_permissions() {
    print_info "Setting permissions..."
    
    # Create user if it doesn't exist
    if ! id "$SERVICE_USER" &>/dev/null; then
        sudo useradd -r -s /bin/false "$SERVICE_USER"
        print_success "Created user: $SERVICE_USER"
    else
        print_info "User already exists: $SERVICE_USER"
    fi
    
    # Set ownership
    sudo chown -R "$SERVICE_USER:$SERVICE_USER" "$DEPLOY_DIR"
    
    # Set permissions
    sudo chmod 755 "$DEPLOY_DIR"
    sudo chmod 644 "$DEPLOY_DIR"/*.jar
    sudo chmod 755 "$LOG_DIR"
    
    print_success "Permissions set"
}

create_systemd_service() {
    print_info "Creating systemd service..."
    
    cat > /tmp/metapilot.service << EOF
[Unit]
Description=MetaPilot Web Application
After=network.target

[Service]
Type=simple
User=$SERVICE_USER
WorkingDirectory=$DEPLOY_DIR
ExecStart=/usr/bin/java -Xms512m -Xmx4g -jar $DEPLOY_DIR/metapilot-web-1.0.0-beta.2.jar
Restart=always
RestartSec=10
StandardOutput=append:$LOG_DIR/stdout.log
StandardError=append:$LOG_DIR/stderr.log

[Install]
WantedBy=multi-user.target
EOF

    sudo cp /tmp/metapilot.service /etc/systemd/system/
    sudo systemctl daemon-reload
    
    print_success "Service created"
}

create_startup_script() {
    print_info "Creating startup script..."
    
    cat > "$DEPLOY_DIR/start.sh" << 'EOF'
#!/bin/bash
# MetaPilot Startup Script

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"

echo "Starting MetaPilot Web Application..."
java -Xms512m -Xmx4g -jar metapilot-web-1.0.0-beta.2.jar

EOF

    sudo chmod +x "$DEPLOY_DIR/start.sh"
    sudo chown "$SERVICE_USER:$SERVICE_USER" "$DEPLOY_DIR/start.sh"
    
    print_success "Startup script created"
}

show_summary() {
    echo ""
    echo -e "${GREEN}================================================${NC}"
    echo -e "${GREEN}  Deployment Complete!${NC}"
    echo -e "${GREEN}================================================${NC}"
    echo ""
    echo "Deployment location: $DEPLOY_DIR"
    echo ""
    echo "To start the service:"
    echo "  sudo systemctl start $SERVICE_NAME"
    echo ""
    echo "To enable on boot:"
    echo "  sudo systemctl enable $SERVICE_NAME"
    echo ""
    echo "To check status:"
    echo "  sudo systemctl status $SERVICE_NAME"
    echo ""
    echo "To view logs:"
    echo "  tail -f $LOG_DIR/stdout.log"
    echo ""
    echo "Or run manually:"
    echo "  sudo -u $SERVICE_USER $DEPLOY_DIR/start.sh"
    echo ""
    echo "The web UI will be available at:"
    echo "  http://localhost:8080"
    echo ""
}

# Main execution
main() {
    print_header
    
    check_jar_exists
    create_directories
    copy_files
    set_permissions
    create_startup_script
    create_systemd_service
    
    show_summary
}

main "$@"

