# MetaPilot HPC Deployment Guide

## Overview

This guide covers deploying MetaPilot on High-Performance Computing (HPC) clusters and virtual machines for production use in research environments.

## Architecture

### Components

MetaPilot HPC consists of:

1. **Web Server**: Spring Boot application with REST API
2. **Workflow Engine**: Executes analysis workflows
3. **Tool Integration**: FragPipe, Sage, DIA-NN
4. **Database**: MAG database management
5. **Results**: Storage and retrieval

### Typical Deployment

```
┌─────────────────────────────────────────┐
│  User's Browser                         │
│  http://vm.internal:8080                │
└─────────────────┬───────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────┐
│  MetaPilot Web Application              │
│  - REST API (Job submission)            │
│  - WebSocket (Progress updates)         │
│  - File management                      │
└─────────────────┬───────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────┐
│  Workflow Execution                     │
│  - FragPipe workflows                   │
│  - Sage workflows                       │
│  - DIA-NN workflows                     │
└─────────────────┬───────────────────────┘
                  │
                  ▼
┌─────────────────────────────────────────┐
│  External Tools                         │
│  /opt/tools/fragpipe                    │
│  /opt/tools/sage                        │
│  /opt/tools/diann                       │
└─────────────────────────────────────────┘
```

## VM Specifications

### Recommended Configuration

For a single-user or small group VM:

- **CPU**: 4-8 cores
- **RAM**: 16-32 GB
- **Storage**: 500 GB+ SSD
- **OS**: Ubuntu 20.04 LTS or CentOS 8

### Multi-User Configuration

For larger groups:

- **CPU**: 8-16 cores
- **RAM**: 32-64 GB
- **Storage**: 1 TB+ SSD
- **Network**: Gigabit connection

## Installation Steps

### 1. Initial VM Setup

```bash
# Update system
sudo apt-get update && sudo apt-get upgrade -y  # Ubuntu
sudo yum update -y  # CentOS

# Install Java 17
sudo apt-get install openjdk-17-jdk maven -y  # Ubuntu
sudo yum install java-17-openjdk-devel maven -y  # CentOS

# Verify
java -version
mvn -version
```

### 2. Install MetaPilot

```bash
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

### 3. Install Search Engines

#### FragPipe

```bash
# Create tools directory
sudo mkdir -p /opt/tools
cd /opt/tools

# Download FragPipe
wget https://github.com/Nesvilab/FragPipe/releases/download/[version]/fragpipe-[version].zip
unzip fragpipe-[version].zip
mv fragpipe-[version] fragpipe

# Set permissions
sudo chmod -R 755 fragpipe
```

#### Sage

```bash
# Download Sage for Linux
wget https://github.com/lazear/sage/releases/download/[version]/sage-linux-x86_64

# Install
sudo mv sage-linux-x86_64 /usr/local/bin/sage
sudo chmod +x /usr/local/bin/sage

# Verify
sage --version
```

#### DIA-NN

```bash
# Download DIA-NN
wget https://github.com/vdemichev/DiaNN/releases/download/[version]/DIA-NN-[version]-Linux.x86_64.zip
unzip DIA-NN-[version]-Linux.x86_64.zip -d /opt/tools/
mv /opt/tools/DIA-NN-[version] /opt/tools/diann

# Make executable
sudo chmod +x /opt/tools/diann/DIA-NN
```

### 4. Configure Resource Paths

Edit `/opt/metapilot/config/meta_resources.params`:

```properties
# Linux Tool Paths
FRAGPIPE_PATH=/opt/tools/fragpipe
SAGE_PATH=/usr/local/bin
DIANN_PATH=/opt/tools/diann

# Shared Directories
MAG_DB_DIR=/opt/metapilot/data/mag_databases
RESULTS_DIR=/opt/metapilot/data/results
TEMP_DIR=/tmp/metapilot

# Python for reports
PYTHON_PATH=/usr/bin/python3
```

### 5. Create Shared Data Directory

```bash
# MAG databases
sudo mkdir -p /opt/metapilot/data/mag_databases
sudo chmod 755 /opt/metapilot/data/mag_databases

# Results
sudo mkdir -p /opt/metapilot/data/results
sudo chmod 755 /opt/metapilot/data/results

# Set ownership
sudo chown -R metapilot:metapilot /opt/metapilot/data
```

### 6. Start Service

```bash
# Start MetaPilot
sudo systemctl start metapilot-web

# Enable on boot
sudo systemctl enable metapilot-web

# Check status
sudo systemctl status metapilot-web
```

## Network Configuration

### Internal Network Access

For institute-only access:

```bash
# Update application.properties
echo "server.address=0.0.0.0" >> metapilot-web/src/main/resources/application.properties

# Configure firewall
sudo ufw allow from 192.168.0.0/16 to any port 8080  # Example internal subnet
```

### Access URL

Users access at: `http://vm-internal-ip:8080`

For domain access, set up DNS:

```bash
# Example: metapilot.internal.yourdomain.edu
# Point to VM IP address
```

## User Workflow

### 1. Access Web UI

Open browser: `http://vm-ip:8080`

### 2. Submit Job

1. Fill in parameter file path
2. Fill in resource file path
3. Click "Submit Job"
4. Receive job ID

### 3. Monitor Progress

- Real-time progress bar
- Status messages
- Job completion notification

### 4. Retrieve Results

Results are saved to: `/opt/metapilot/data/results/{job-id}/`

## MAG Database Management

### Download MAG Databases

Through web UI or CLI:

1. Access MAG management page
2. Browse available databases
3. Select and download
4. Verify installation

### Shared Database Location

All MAG databases stored at:
```
/opt/metapilot/data/mag_databases/
├── mgyg-catalogue-v1/
├── mgyg-catalogue-v2/
└── ...
```

## Performance Tuning

### JVM Options

For large workloads:

```bash
# Edit systemd service
sudo nano /etc/systemd/system/metapilot.service

# Update ExecStart
ExecStart=/usr/bin/java -Xms2g -Xmx16g -XX:+UseG1GC -jar ...
```

### Concurrent Jobs

Limit in `WorkflowService.java`:

```java
executorService = Executors.newFixedThreadPool(4); // Adjust based on CPU cores
```

### Disk I/O

Use SSD for:
- `/tmp/metapilot`
- `/opt/metapilot/data/results`
- Tool working directories

## Monitoring

### Log Files

```bash
# Application logs
tail -f /opt/metapilot/logs/stdout.log

# Errors
tail -f /opt/metapilot/logs/stderr.log

# System logs
journalctl -u metapilot-web -f
```

### Resource Monitoring

```bash
# CPU/Memory usage
htop

# Disk usage
df -h /opt/metapilot

# Running processes
ps aux | grep java
```

### Health Checks

```bash
# API status
curl http://localhost:8080/api/jobs/status

# Check process
ps aux | grep metapilot
```

## Backup Strategy

### Configuration Backup

```bash
# Daily backup
sudo tar -czf /backup/metapilot-config-$(date +%Y%m%d).tar.gz /opt/metapilot/config
```

### Results Backup

```bash
# Weekly backup
sudo rsync -av /opt/metapilot/data/results/ /backup/results/
```

### Database Backup

MAG databases can be large. Store externally or compress:

```bash
# Compress MAG databases
sudo tar -czf /backup/mag-db-$(date +%Y%m%d).tar.gz /opt/metapilot/data/mag_databases
```

## Maintenance

### Regular Updates

```bash
# Pull latest code
cd /path/to/metapilot
git pull origin main

# Rebuild
./build_linux.sh

# Restart service
sudo systemctl restart metapilot-web
```

### Log Rotation

Configure logrotate:

```bash
sudo nano /etc/logrotate.d/metapilot

# Add:
/opt/metapilot/logs/*.log {
    daily
    rotate 7
    compress
    delaycompress
    missingok
    notifempty
    create 0640 metapilot metapilot
}
```

### Cleanup Old Results

```bash
# Automated cleanup script
find /opt/metapilot/data/results -type d -mtime +90 -exec rm -rf {} \;
```

## Security

### User Access

Currently, MetaPilot has no authentication. For production:

1. Add user authentication (planned)
2. Network firewall rules
3. IP whitelisting
4. HTTPS via reverse proxy

### Firewall Rules

```bash
# Restrict to internal network
sudo ufw deny from any to any port 8080
sudo ufw allow from 192.168.1.0/24 to any port 8080
sudo ufw reload
```

### SSL/HTTPS (via Nginx)

Set up reverse proxy with SSL:

```nginx
server {
    listen 443 ssl http2;
    server_name metapilot.yourdomain.edu;
    
    ssl_certificate /etc/letsencrypt/live/metapilot.yourdomain.edu/fullchain.pem;
    ssl_certificate_key /etc/letsencrypt/live/metapilot.yourdomain.edu/privkey.pem;
    
    location / {
        proxy_pass http://localhost:8080;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection 'upgrade';
        proxy_set_header Host $host;
    }
}
```

## Troubleshooting

### Service Won't Start

```bash
# Check logs
journalctl -u metapilot-web -n 50

# Check Java
java -version

# Check ports
sudo netstat -tulpn | grep 8080
```

### Out of Memory

Increase heap in systemd service:

```bash
sudo nano /etc/systemd/system/metapilot.service
# Change: -Xmx8g to -Xmx16g
sudo systemctl daemon-reload
sudo systemctl restart metapilot-web
```

### WebSocket Issues

Check browser console for errors. Verify:
- Port 8080 is accessible
- Firewall allows connections
- Proxy settings if behind proxy

### Tool Execution Fails

```bash
# Check tool paths in config
cat /opt/metapilot/config/meta_resources.params

# Verify tools are executable
ls -l /opt/tools/fragpipe
ls -l /usr/local/bin/sage
ls -l /opt/tools/diann/DIA-NN

# Test manually
/opt/tools/fragpipe --version
sage --version
/opt/tools/diann/DIA-NN
```

## Scalability

### Horizontal Scaling

For multiple VMs:

1. Deploy MetaPilot on each VM
2. Use load balancer (Nginx/HAProxy)
3. Share results via network storage (NFS)
4. Centralized logging

### Vertical Scaling

For single powerful VM:

1. Increase JVM heap size
2. Add more CPU cores
3. Use faster storage (NVMe SSD)
4. Increase concurrent job pool

## Support

### Getting Help

- Check logs: `/opt/metapilot/logs/`
- Review documentation: `docs/`
- Test with sample data
- Contact: MetaPilot development team

### Common Issues

- **Port conflicts**: Change port in `application.properties`
- **Permission errors**: Check ownership with `ls -la`
- **Tool not found**: Verify paths in config file
- **Memory errors**: Increase `-Xmx` parameter

## Appendix

### Quick Commands

```bash
# Start
sudo systemctl start metapilot-web

# Stop
sudo systemctl stop metapilot-web

# Restart
sudo systemctl restart metapilot-web

# Status
sudo systemctl status metapilot-web

# Logs
tail -f /opt/metapilot/logs/stdout.log

# Configuration
nano /opt/metapilot/config/meta_resources.params

# Rebuild
cd /path/to/metapilot && ./build_linux.sh

# Deploy
sudo ./deploy_linux.sh
```

### Directory Structure

```
/opt/metapilot/
├── metapilot-web-1.0.0-beta.2.jar
├── start.sh
├── config/
│   └── meta_resources.params
├── logs/
│   ├── stdout.log
│   └── stderr.log
└── data/
    ├── mag_databases/
    └── results/
```

---

**Version**: 1.0.0-beta.3  
**Platform**: Linux HPC/VMs  
**Last Updated**: 2024

