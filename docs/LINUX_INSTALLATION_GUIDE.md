# MetaPilot Linux Installation Guide

## Overview

This guide provides step-by-step instructions for installing and running MetaPilot on Linux systems, including HPC environments and virtual machines.

## Prerequisites

### System Requirements

- **Operating System**: Linux (Ubuntu 18.04+, CentOS 7+, RHEL 8+, Debian 10+)
- **Java**: JDK 17 or higher
- **Memory**: Minimum 4GB RAM (8GB+ recommended)
- **Disk Space**: Minimum 10GB free space
- **Network**: Internet connection for initial setup

### Software Dependencies

- **Java Development Kit (JDK)**: Required for running MetaPilot
- **Maven**: Required for building from source
- **Git**: Required for cloning repository (if installing from source)

## Installation Methods

### Method 1: Pre-built JAR (Recommended)

If you have a pre-built JAR file:

```bash
# Download or copy the JAR file
wget https://example.com/metapilot-web-1.0.0-beta.2.jar

# Create installation directory
sudo mkdir -p /opt/metapilot
sudo cp metapilot-web-1.0.0-beta.2.jar /opt/metapilot/

# Set permissions
sudo useradd -r -s /bin/false metapilot
sudo chown metapilot:metapilot /opt/metapilot
```

### Method 2: Build from Source

If building from source:

```bash
# Install prerequisites
sudo apt-get update  # Ubuntu/Debian
sudo yum update      # CentOS/RHEL

# Install Java 17
sudo apt-get install openjdk-17-jdk  # Ubuntu/Debian
sudo yum install java-17-openjdk-devel  # CentOS/RHEL

# Install Maven
sudo apt-get install maven  # Ubuntu/Debian
sudo yum install maven  # CentOS/RHEL

# Clone repository
git clone https://github.com/your-org/metapilot.git
cd metapilot

# Build
chmod +x build_linux.sh
./build_linux.sh

# Deploy
chmod +x deploy_linux.sh
sudo ./deploy_linux.sh
```

## Configuration

### 1. Create Configuration Directory

```bash
sudo mkdir -p /opt/metapilot/config
sudo mkdir -p /opt/metapilot/logs
sudo mkdir -p /opt/metapilot/data
sudo chown -R metapilot:metapilot /opt/metapilot
```

### 2. Configure Resources

Create or edit `/opt/metapilot/config/meta_resources.params`:

```properties
# Tool Paths for Linux
# Update these paths to match your installation

# FragPipe (Linux)
FRAGPIPE_PATH=/opt/tools/fragpipe

# Sage (Linux)
SAGE_PATH=/opt/tools/sage

# DIA-NN (Linux)
DIANN_PATH=/opt/tools/diann

# Python (for report generation)
PYTHON_PATH=/usr/bin/python3

# Temporary directory
TEMP_DIR=/tmp/metapilot
```

### 3. Environment Variables

Optionally set environment variables:

```bash
# Add to ~/.bashrc or /etc/profile
export METAPILOT_HOME=/opt/metapilot
export METAPILOT_CONFIG=/opt/metapilot/config
export JAVA_OPTS="-Xms512m -Xmx4g"
```

## Running MetaPilot

### Option 1: Manual Start

```bash
cd /opt/metapilot
java -Xms512m -Xmx4g -jar metapilot-web-1.0.0-beta.2.jar
```

Access the web UI at: `http://localhost:8080`

### Option 2: Systemd Service

If you used `deploy_linux.sh`, the systemd service is already configured:

```bash
# Start service
sudo systemctl start metapilot-web

# Enable on boot
sudo systemctl enable metapilot-web

# Check status
sudo systemctl status metapilot-web

# View logs
tail -f /opt/metapilot/logs/stdout.log
```

### Option 3: Screen/Tmux Session

For HPC environments without systemd:

```bash
# Start with screen
screen -S metapilot
cd /opt/metapilot
java -Xms512m -Xmx4g -jar metapilot-web-1.0.0-beta.2.jar

# Detach: Ctrl+A, then D
# Reattach: screen -r metapilot
```

Or with tmux:

```bash
tmux new -s metapilot
java -Xms512m -Xmx4g -jar metapilot-web-1.0.0-beta.2.jar

# Detach: Ctrl+B, then D
# Reattach: tmux attach -t metapilot
```

## Installing Search Engines

### FragPipe (Linux)

FragPipe has native Linux support:

```bash
# Download FragPipe
wget https://github.com/Nesvilab/FragPipe/releases/download/[version]/fragpipe-[version].zip

# Extract
unzip fragpipe-[version].zip -d /opt/tools/
cd /opt/tools/fragpipe

# Make executable
chmod +x fragpipe

# Test
./fragpipe --version
```

### Sage (Linux)

Sage has native Linux support:

```bash
# Download Sage
wget https://github.com/lazear/sage/releases/download/[version]/sage-linux-x86_64

# Install
sudo mv sage-linux-x86_64 /usr/local/bin/sage
sudo chmod +x /usr/local/bin/sage

# Test
sage --version
```

### DIA-NN (Linux)

DIA-NN has native Linux support:

```bash
# Download DIA-NN
wget https://github.com/vdemichev/DiaNN/releases/download/[version]/DIA-NN-[version]-Linux.x86_64.zip

# Extract
unzip DIA-NN-[version]-Linux.x86_64.zip -d /opt/tools/
cd /opt/tools/diann

# Make executable
chmod +x DIA-NN

# Test
./DIA-NN
```

## Port Configuration

By default, MetaPilot runs on port 8080. To change this:

Edit `metapilot-web/src/main/resources/application.properties`:

```properties
server.port=8080  # Change to desired port
server.address=0.0.0.0  # Listen on all interfaces
```

Or set via environment variable:

```bash
export SERVER_PORT=8080
java -jar metapilot-web-1.0.0-beta.2.jar
```

Or via command line:

```bash
java -jar metapilot-web-1.0.0-beta.2.jar --server.port=8080
```

## Firewall Configuration

For HPC environments, you may need to configure firewall rules:

```bash
# Ubuntu/Debian (UFW)
sudo ufw allow 8080/tcp

# CentOS/RHEL (firewalld)
sudo firewall-cmd --permanent --add-port=8080/tcp
sudo firewall-cmd --reload

# iptables (generic)
sudo iptables -A INPUT -p tcp --dport 8080 -j ACCEPT
```

## Verification

### 1. Check Java Version

```bash
java -version
# Should show: openjdk version "17.x.x" or higher
```

### 2. Check Application Status

```bash
# Via systemd
sudo systemctl status metapilot-web

# Via curl
curl http://localhost:8080/api/jobs/status

# Expected response: {"status":"running"}
```

### 3. Access Web UI

Open browser: `http://your-server-ip:8080`

You should see the MetaPilot web interface.

## Troubleshooting

### Port Already in Use

**Symptom**: `Address already in use`

**Solution**:
```bash
# Find process using port 8080
sudo lsof -i :8080
sudo netstat -tulpn | grep 8080

# Kill process or change MetaPilot port
```

### Java Not Found

**Symptom**: `java: command not found`

**Solution**:
```bash
# Install Java
sudo apt-get install openjdk-17-jdk  # Ubuntu/Debian
sudo yum install java-17-openjdk-devel  # CentOS/RHEL

# Set JAVA_HOME
export JAVA_HOME=/usr/lib/jvm/java-17-openjdk
```

### Permission Denied

**Symptom**: `Permission denied` when starting

**Solution**:
```bash
# Set correct permissions
sudo chown -R metapilot:metapilot /opt/metapilot
sudo chmod +x metapilot-web-1.0.0-beta.2.jar
```

### Out of Memory

**Symptom**: `OutOfMemoryError`

**Solution**:
```bash
# Increase heap size
java -Xms1g -Xmx8g -jar metapilot-web-1.0.0-beta.2.jar
```

### WebSocket Connection Failed

**Symptom**: WebSocket not connecting in browser

**Solution**:
- Check firewall allows port 8080
- Verify SockJS CDN is accessible
- Check browser console for errors

## Security Considerations

### Production Deployment

For production use:

1. **Change Default Port**: Don't use 8080 on public servers
2. **Enable HTTPS**: Use reverse proxy (Nginx/Apache) with SSL
3. **Authentication**: Add user authentication
4. **Firewall**: Only allow necessary IPs
5. **Updates**: Keep Java and MetaPilot updated
6. **Monitoring**: Set up log rotation and monitoring

### Reverse Proxy with Nginx

Example Nginx configuration:

```nginx
server {
    listen 80;
    server_name metapilot.your-domain.com;
    
    location / {
        proxy_pass http://localhost:8080;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection 'upgrade';
        proxy_set_header Host $host;
        proxy_cache_bypass $http_upgrade;
    }
}
```

## Next Steps

After installation:

1. Read the [HPC Deployment Guide](HPC_DEPLOYMENT_GUIDE.md)
2. Review the [CLI Reference Manual](CLI_REFERENCE_MANUAL.md)
3. Configure workflow tools (FragPipe, Sage, DIA-NN)
4. Test with sample data
5. Set up monitoring and backups

## Support

For issues or questions:

- Check logs: `/opt/metapilot/logs/`
- Review documentation: `docs/`
- Contact: MetaPilot development team

## Appendix

### Quick Reference

```bash
# Build
./build_linux.sh

# Deploy
sudo ./deploy_linux.sh

# Start
sudo systemctl start metapilot-web

# Stop
sudo systemctl stop metapilot-web

# Restart
sudo systemctl restart metapilot-web

# Logs
tail -f /opt/metapilot/logs/stdout.log

# Status
curl http://localhost:8080/api/jobs/status
```

---

**Version**: 1.0.0-beta.3  
**Last Updated**: 2024  
**Platform**: Linux (x86_64)

