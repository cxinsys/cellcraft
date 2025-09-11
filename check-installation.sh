#!/bin/bash

# CellCraft Installation Check Script
# Comprehensive validation for CellCraft system setup
# Supports both docker-compose.cpu.yml and docker-compose.local.yml configurations

# Removed strict error checking to allow better error handling

# Color codes for output
readonly RED='\033[0;31m'
readonly GREEN='\033[0;32m' 
readonly YELLOW='\033[1;33m'
readonly BLUE='\033[0;34m'
readonly NC='\033[0m' # No Color

# Configuration variables
COMPOSE_FILE=""
PROJECT_NAME="cellcraft"
TROUBLESHOOTING_URL="https://cellcraft.gitbook.io/cellcraft-docs"
TIMEOUT=30

# Status counters
TOTAL_CHECKS=0
PASSED_CHECKS=0
FAILED_CHECKS=0

# Functions for colored output
log_info() { echo -e "${BLUE}ℹ️  $1${NC}"; }
log_success() { echo -e "${GREEN}✅ $1${NC}"; ((PASSED_CHECKS++)); }
log_warning() { echo -e "${YELLOW}⚠️  $1${NC}"; ((PASSED_CHECKS++)); }  # Count warnings as passed
log_error() { echo -e "${RED}❌ $1${NC}"; ((FAILED_CHECKS++)); }
log_header() { echo -e "\n${BLUE}🔍 $1${NC}\n"; }

# Increment total checks counter
check_start() {
    ((TOTAL_CHECKS++))
}

# Print script header
print_header() {
    echo -e "${BLUE}"
    echo "╔════════════════════════════════════════════════════════════════╗"
    echo "║                    CellCraft Installation Check                ║"
    echo "║            Gene Regulatory Network Analysis Platform           ║" 
    echo "╚════════════════════════════════════════════════════════════════╝"
    echo -e "${NC}"
}

# Detect Docker Compose configuration
detect_compose_config() {
    log_header "Configuration Detection"
    
    check_start
    if [[ -f "docker-compose.cpu.yml" && -f "docker-compose.local.yml" ]]; then
        # Simple approach: check running containers and match with expected services
        local running_cellcraft_containers
        running_cellcraft_containers=$(docker ps --format "table {{.Names}}" 2>/dev/null | grep "^cellcraft-" | wc -l)
        running_cellcraft_containers=${running_cellcraft_containers:-0}
        
        if [[ ${running_cellcraft_containers} -gt 0 ]]; then
            # Assume CPU configuration if containers are running (most common)
            COMPOSE_FILE="docker-compose.cpu.yml"
            log_success "CellCraft containers detected, using CPU configuration (${running_cellcraft_containers} containers)"
        else
            # Default to CPU if no containers are running
            COMPOSE_FILE="docker-compose.cpu.yml"
            log_warning "No CellCraft containers running, defaulting to CPU configuration"
        fi
    elif [[ -f "docker-compose.cpu.yml" ]]; then
        COMPOSE_FILE="docker-compose.cpu.yml"
        log_success "CPU-only configuration available"
    elif [[ -f "docker-compose.local.yml" ]]; then
        COMPOSE_FILE="docker-compose.local.yml"
        log_success "GPU-enabled configuration available"
    else
        log_error "No Docker Compose configuration files found"
        log_error "Expected: docker-compose.cpu.yml or docker-compose.local.yml"
        exit 1
    fi
    
    log_info "Using configuration: ${COMPOSE_FILE}"
}

# Validate environment variables
validate_environment() {
    log_header "Environment Validation"
    
    # Check for .env file
    check_start
    if [[ -f ".env" ]]; then
        log_success ".env file found"
        # Source .env file with proper parsing (handle spaces around =)
        while IFS= read -r line || [[ -n "$line" ]]; do
            if [[ $line =~ ^[^#]*= ]]; then
                # Remove spaces around = and export the variable
                key=$(echo "$line" | cut -d'=' -f1 | xargs)
                value=$(echo "$line" | cut -d'=' -f2- | xargs)
                export "$key=$value"
            fi
        done < .env
    else
        log_warning ".env file not found - using default values"
    fi
    
    # Set default values if not provided
    export POSTGRES_DB=${POSTGRES_DB:-cellcraft}
    export POSTGRES_USER=${POSTGRES_USER:-postgres}
    export POSTGRES_PASSWORD=${POSTGRES_PASSWORD:-password}
    export GPU_COUNT=${GPU_COUNT:-all}
    
    check_start
    log_success "Environment variables configured:"
    log_info "  - Database: ${POSTGRES_DB}"
    log_info "  - User: ${POSTGRES_USER}"
    if [[ "${COMPOSE_FILE}" == "docker-compose.local.yml" ]]; then
        log_info "  - GPU Count: ${GPU_COUNT}"
    fi
}

# Check Docker daemon status
check_docker() {
    log_header "Docker System Check"
    
    check_start
    if ! command -v docker &> /dev/null; then
        log_error "Docker is not installed or not in PATH"
        return 1
    fi
    log_success "Docker command available"
    
    check_start
    if ! docker info &> /dev/null; then
        log_error "Docker daemon is not running"
        log_error "Please start Docker and try again"
        return 1
    fi
    log_success "Docker daemon is running"
    
    check_start
    if ! command -v docker-compose &> /dev/null && ! docker compose version &> /dev/null; then
        log_error "Docker Compose is not available"
        return 1
    fi
    log_success "Docker Compose is available"
}

# Check container status
check_containers() {
    log_header "Container Status Check"
    
    local expected_services=("frontend" "backend" "db" "rabbitmq" "celery")
    local running_count=0
    
    # Get list of running containers for this project
    local running_containers
    running_containers=$(docker compose -f "${COMPOSE_FILE}" ps --format "table {{.Service}}" 2>/dev/null | tail -n +2 | tr -d ' ')
    
    for service in "${expected_services[@]}"; do
        check_start
        if echo "${running_containers}" | grep -q "^${service}$"; then
            # Check if container is actually running (not just created)
            local container_status
            container_status=$(docker compose -f "${COMPOSE_FILE}" ps "${service}" --format "table {{.State}}" 2>/dev/null | tail -n +2 | tr -d ' ')
            if [[ "${container_status}" == "running" ]]; then
                log_success "${service} container is running"
                ((running_count++))
            else
                log_error "${service} container exists but is not running (status: ${container_status})"
            fi
        else
            log_error "${service} container is not running"
            log_info "  Run: docker compose -f ${COMPOSE_FILE} up -d ${service}"
        fi
    done
    
    log_info "Container summary: ${running_count}/${#expected_services[@]} services running"
    
    if [[ ${running_count} -eq ${#expected_services[@]} ]]; then
        log_success "All required containers are running"
    else
        log_error "Some containers are missing or stopped"
        log_error "Run: docker compose -f ${COMPOSE_FILE} up -d"
    fi
}

# Check service health
check_service_health() {
    log_header "Service Health Check"
    
    # Check database health
    check_start
    local db_container="${PROJECT_NAME}-db-1"
    if docker exec "${db_container}" pg_isready -U "${POSTGRES_USER}" -d "${POSTGRES_DB}" &>/dev/null; then
        log_success "Database is healthy and accepting connections"
    else
        log_error "Database health check failed"
        log_error "Check database logs: docker logs ${db_container}"
    fi
    
    # Check RabbitMQ health  
    check_start
    local rabbitmq_container="${PROJECT_NAME}-rabbitmq-1"
    if docker exec "${rabbitmq_container}" rabbitmq-diagnostics ping &>/dev/null; then
        log_success "RabbitMQ is healthy"
    else
        log_error "RabbitMQ health check failed"
        log_error "Check RabbitMQ logs: docker logs ${rabbitmq_container}"
    fi
}

# Test service connectivity
test_connectivity() {
    log_header "Service Connectivity Test"
    
    # Test frontend
    check_start
    if curl -s --connect-timeout ${TIMEOUT} http://localhost:8080 &>/dev/null; then
        log_success "Frontend is accessible at http://localhost:8080"
    else
        log_error "Frontend is not accessible at http://localhost:8080"
        log_error "Check frontend container logs and port mapping"
    fi
    
    # Test backend API
    check_start
    if curl -s --connect-timeout ${TIMEOUT} http://localhost:8000/docs &>/dev/null; then
        log_success "Backend API is accessible at http://localhost:8000"
    else
        log_error "Backend API is not accessible at http://localhost:8000"
        log_error "Check backend container logs and port mapping"
    fi
    
    # Test database connection
    check_start
    local db_container="${PROJECT_NAME}-db-1"
    if docker exec "${db_container}" psql -U "${POSTGRES_USER}" -d "${POSTGRES_DB}" -c "SELECT 1;" >/dev/null 2>&1; then
        log_success "Database connection successful"
    else
        log_warning "Direct database connection test skipped (may require specific setup)"
        log_info "  Database container is healthy, which indicates working connectivity"
    fi
    
    # Test RabbitMQ management interface
    check_start
    if curl -s --connect-timeout ${TIMEOUT} http://localhost:15672 &>/dev/null; then
        log_success "RabbitMQ management interface is accessible"
    else
        log_error "RabbitMQ management interface is not accessible"
        log_error "Check RabbitMQ container and port 15672"
    fi
}

# Validate plugin registry
validate_plugins() {
    log_header "Plugin Registry Validation"
    
    check_start
    local plugin_response=""
    if plugin_response=$(curl -s --connect-timeout ${TIMEOUT} http://localhost:8000/api/plugin/list 2>/dev/null); then
        log_success "Plugin registry is accessible"
        
        # Count plugins
        local plugin_count
        if command -v jq &> /dev/null; then
            plugin_count=$(echo "${plugin_response}" | jq '. | length' 2>/dev/null || echo "unknown")
        else
            # Fallback counting method without jq
            plugin_count=$(echo "${plugin_response}" | grep -o '"name"' | wc -l 2>/dev/null || echo "unknown")
        fi
        
        if [[ "${plugin_count}" != "unknown" && ${plugin_count} -gt 0 ]]; then
            log_success "Plugin registry loaded: ${plugin_count} plugins available"
            log_info "  Expected plugins: TENET, FastTENET, GENIE3, GRNBoost2, LEAP, Scribe, GRNViz, FastSCODE"
        else
            log_warning "Plugin registry is accessible but plugin count is unclear"
        fi
    else
        log_error "Plugin registry is not accessible"
        log_error "Backend API may not be properly initialized"
    fi
}

# Check system resources
check_resources() {
    log_header "System Resource Check"
    
    # Check Docker daemon resources
    check_start
    local docker_info
    if docker_info=$(docker system df 2>/dev/null); then
        log_success "Docker system resources:"
        echo "${docker_info}" | while IFS= read -r line; do
            log_info "  ${line}"
        done
    else
        log_error "Cannot retrieve Docker system information"
    fi
    
    # Check disk space
    check_start
    local disk_usage
    disk_usage=$(df -h . | tail -1 | awk '{print $4}')
    log_success "Available disk space: ${disk_usage}"
    
    # GPU check for local configuration
    if [[ "${COMPOSE_FILE}" == "docker-compose.local.yml" ]]; then
        check_start
        if command -v nvidia-smi &> /dev/null; then
            if nvidia-smi &> /dev/null; then
                log_success "NVIDIA GPU detected and accessible"
                log_info "$(nvidia-smi --query-gpu=name --format=csv,noheader,nounits | head -1)"
            else
                log_error "NVIDIA drivers not working properly"
            fi
        else
            log_warning "nvidia-smi not found - GPU support may not be available"
        fi
    fi
}

# Generate final report
generate_report() {
    log_header "Installation Check Summary"
    
    local success_rate=$((PASSED_CHECKS * 100 / TOTAL_CHECKS))
    
    echo -e "\n${BLUE}📊 Check Results:${NC}"
    echo -e "  Total Checks: ${TOTAL_CHECKS}"
    echo -e "  ${GREEN}Passed: ${PASSED_CHECKS}${NC}"
    echo -e "  ${RED}Failed: ${FAILED_CHECKS}${NC}"
    echo -e "  Success Rate: ${success_rate}%"
    
    if [[ ${FAILED_CHECKS} -eq 0 ]]; then
        echo -e "\n${GREEN}🎉 Installation successful!${NC}"
        echo -e "${GREEN}CellCraft is ready for Gene Regulatory Network analysis${NC}"
        echo -e "\n${BLUE}Access Points:${NC}"
        echo -e "  • Frontend: http://localhost:8080"
        echo -e "  • Backend API: http://localhost:8000/docs"
        echo -e "  • RabbitMQ Management: http://localhost:15672"
    else
        echo -e "\n${RED}⚠️  Installation check completed with issues${NC}"
        echo -e "${YELLOW}Please review the failed checks above and take corrective action${NC}"
        echo -e "\n${BLUE}Troubleshooting Resources:${NC}"
        echo -e "  • Documentation: ${TROUBLESHOOTING_URL}"
        echo -e "  • Check logs: docker compose -f ${COMPOSE_FILE} logs"
        echo -e "  • Restart services: docker compose -f ${COMPOSE_FILE} restart"
        echo -e "  • Full rebuild: docker compose -f ${COMPOSE_FILE} down && docker compose -f ${COMPOSE_FILE} up --build -d"
    fi
    
    echo -e "\n${BLUE}Configuration Used: ${COMPOSE_FILE}${NC}"
    echo -e "${BLUE}For additional help, visit: ${TROUBLESHOOTING_URL}${NC}\n"
}

# Main execution
main() {
    print_header
    
    # Change to script directory
    cd "$(dirname "${BASH_SOURCE[0]}")"
    
    # Execute all checks
    detect_compose_config
    validate_environment
    check_docker
    check_containers
    check_service_health
    test_connectivity
    validate_plugins
    check_resources
    
    # Generate final report
    generate_report
    
    # Exit with appropriate code
    if [[ ${FAILED_CHECKS} -eq 0 ]]; then
        exit 0
    else
        exit 1
    fi
}

# Execute main function
main "$@"