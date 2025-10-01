#!/usr/bin/env julia

"""
Clean Code Analysis Script for SparseQEEcNEO.jl

This script analyzes the codebase for Clean Code compliance:
- Functions longer than 20 lines
- Complex cyclomatic complexity 
- Poor naming conventions
- Magic numbers
- Deep nesting levels
"""

using Pkg

function count_lines_in_function(lines, start_idx, func_name)
    """Count lines in a function definition"""
    line_count = 0
    brace_count = 0
    in_function = false
    
    for i in start_idx:length(lines)
        line = strip(lines[i])
        
        # Skip empty lines and comments in count
        if isempty(line) || startswith(line, "#")
            continue
        end
        
        # Track function boundaries
        if occursin("function", line) || occursin("$(func_name)(", line)
            in_function = true
            brace_count += count(c -> c == '{', line) - count(c -> c == '}', line)
            line_count += 1
        elseif in_function
            line_count += 1
            
            # Track braces/keywords for function end
            if occursin("end", line) && !occursin("endof", line) && !occursin("endswith", line)
                break
            end
            
            # Alternative: look for return or function end patterns
            if occursin("return", line) && i < length(lines) - 1
                next_line = strip(lines[i + 1])
                if isempty(next_line) || startswith(next_line, "function") || startswith(next_line, "end")
                    break
                end
            end
        end
    end
    
    return line_count
end

function analyze_function_complexity(lines, start_idx, func_name)
    """Analyze cyclomatic complexity of a function"""
    complexity = 1  # Base complexity
    nesting_level = 0
    max_nesting = 0
    
    for i in start_idx:length(lines)
        line = strip(lines[i])
        
        if isempty(line) || startswith(line, "#")
            continue
        end
        
        # Count complexity-adding constructs
        complexity += count(x -> occursin(x, line), ["if ", "elseif", "else", "for ", "while ", "try", "catch"])
        complexity += count(x -> occursin(x, line), ["&&", "||"])
        
        # Track nesting level
        if occursin("if ", line) || occursin("for ", line) || occursin("while ", line) || occursin("try", line)
            nesting_level += 1
            max_nesting = max(max_nesting, nesting_level)
        elseif occursin("end", line)
            nesting_level = max(0, nesting_level - 1)
        end
        
        # Break at function end
        if occursin("end", line) && nesting_level == 0
            break
        end
    end
    
    return complexity, max_nesting
end

function find_magic_numbers(content)
    """Find potential magic numbers in code"""
    magic_numbers = Set{String}()
    
    # Look for numeric literals (excluding common ones like 0, 1, 2)
    number_pattern = r"\b(\d+\.?\d*)\b"
    
    for match in eachmatch(number_pattern, content)
        num_str = match.captures[1]
        num = tryparse(Float64, num_str)
        
        if num !== nothing && !(num in [0, 1, 2, 3, 10, 100, 1000])
            push!(magic_numbers, num_str)
        end
    end
    
    return collect(magic_numbers)
end

function analyze_file(filepath)
    """Analyze a single Julia file for Clean Code violations"""
    println("\n" * "="^60)
    println("Analyzing: $(basename(filepath))")
    println("="^60)
    
    content = read(filepath, String)
    lines = split(content, '\n')
    
    violations = []
    function_pattern = r"function\s+([a-zA-Z_][a-zA-Z0-9_!]*)"
    
    # Find all functions
    for (i, line) in enumerate(lines)
        line = strip(line)
        
        # Skip comments and empty lines
        if isempty(line) || startswith(line, "#")
            continue
        end
        
        # Look for function definitions
        func_match = match(function_pattern, line)
        if func_match !== nothing
            func_name = func_match.captures[1]
            
            # Count function lines
            func_lines = count_lines_in_function(lines, i, func_name)
            
            # Analyze complexity
            complexity, max_nesting = analyze_function_complexity(lines, i, func_name)
            
            # Check for violations
            issues = String[]
            
            if func_lines > 20
                push!(issues, "Long function ($(func_lines) lines)")
            end
            
            if complexity > 10
                push!(issues, "High complexity ($(complexity))")
            end
            
            if max_nesting > 4
                push!(issues, "Deep nesting ($(max_nesting) levels)")
            end
            
            # Check naming conventions
            if !islowercase(func_name[1]) || occursin("_", func_name)
                if !(func_name in ["run_cneo_hf", "run_cneo_mp2", "sparse_qee_cneo", "sparse_qee_with_cneo"])
                    push!(issues, "Poor naming convention")
                end
            end
            
            if !isempty(issues)
                push!(violations, (func_name, func_lines, complexity, max_nesting, issues))
                
                println("❌ Function: $func_name")
                println("   Lines: $func_lines | Complexity: $complexity | Max nesting: $max_nesting")
                println("   Issues: $(join(issues, ", "))")
            else
                println("✅ Function: $func_name ($(func_lines) lines, complexity: $(complexity))")
            end
        end
    end
    
    # Find magic numbers
    magic_nums = find_magic_numbers(content)
    if !isempty(magic_nums)
        println("\n🔢 Magic numbers found: $(join(magic_nums, ", "))")
    end
    
    return violations
end

function main()
    println("🧹 Clean Code Analysis for SparseQEEcNEO.jl")
    println("="^60)
    
    src_dir = joinpath(pwd(), "src")
    julia_files = [joinpath(src_dir, f) for f in readdir(src_dir) if endswith(f, ".jl")]
    
    total_violations = 0
    all_violations = []
    
    for file in julia_files
        violations = analyze_file(file)
        total_violations += length(violations)
        append!(all_violations, violations)
    end
    
    println("\n" * "="^60)
    println("📊 SUMMARY")
    println("="^60)
    println("Total violations: $total_violations")
    println("Files analyzed: $(length(julia_files))")
    
    if total_violations == 0
        println("🎉 100% Clean Code compliant!")
    else
        println("\n📝 Violations by category:")
        
        long_functions = [v for v in all_violations if any(occursin("Long function", issue) for issue in v[5])]
        complex_functions = [v for v in all_violations if any(occursin("High complexity", issue) for issue in v[5])]
        deep_nesting = [v for v in all_violations if any(occursin("Deep nesting", issue) for issue in v[5])]
        naming_issues = [v for v in all_violations if any(occursin("Poor naming", issue) for issue in v[5])]
        
        println("   Long functions (>20 lines): $(length(long_functions))")
        println("   High complexity (>10): $(length(complex_functions))")
        println("   Deep nesting (>4 levels): $(length(deep_nesting))")
        println("   Naming issues: $(length(naming_issues))")
        
        println("\n🎯 Priority functions to refactor:")
        sorted_violations = sort(all_violations, by = x -> x[2], rev = true)  # Sort by line count
        for (i, (name, lines, complexity, nesting, issues)) in enumerate(sorted_violations[1:min(10, length(sorted_violations))])
            println("   $i. $name ($(lines) lines, complexity: $complexity)")
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end