# Top-level: only find Git once
find_program(GIT_EXECUTABLE git)

# Define reusable function
function(get_git_info dir prefix)
    set(${prefix}_BRANCH "unknown" PARENT_SCOPE)
    set(${prefix}_REVISION "unknown" PARENT_SCOPE)

    if(GIT_EXECUTABLE)
        execute_process(
                COMMAND ${GIT_EXECUTABLE} rev-parse --is-inside-work-tree
                WORKING_DIRECTORY ${dir}
                OUTPUT_VARIABLE _git_check_output     # suppress "true"
                RESULT_VARIABLE IS_GIT_REPO
                OUTPUT_STRIP_TRAILING_WHITESPACE
        )

        if(IS_GIT_REPO EQUAL 0)
            execute_process(
                    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
                    WORKING_DIRECTORY ${dir}
                    OUTPUT_VARIABLE GIT_BRANCH
                    OUTPUT_STRIP_TRAILING_WHITESPACE
            )

            execute_process(
                    COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
                    WORKING_DIRECTORY ${dir}
                    OUTPUT_VARIABLE GIT_COMMIT
                    OUTPUT_STRIP_TRAILING_WHITESPACE
            )

            set(${prefix}_BRANCH "${GIT_BRANCH}" PARENT_SCOPE)
            set(${prefix}_REVISION "${GIT_COMMIT}" PARENT_SCOPE)
        else()
            message(WARNING "${dir} is not a Git repository")
        endif()
    else()
        message(WARNING "Git executable not found")
    endif()
endfunction()


function(ListPluginsGit)
    set(MAX_NAME_WIDTH 25)
    set(MAX_BRANCH_WIDTH 30)
    get_property(plugins GLOBAL PROPERTY ugPluginNames)
    foreach(plugin ${plugins})
        if(${${plugin}} STREQUAL "ON")
            #set(enabledPluginsStr ${enabledPluginsStr} ${plugin})
            get_git_info(${CMAKE_SOURCE_DIR}/plugins/${plugin} "${plugin}_GIT")
            pad_right(${plugin} ${MAX_NAME_WIDTH} plugin_padded)
            pad_right(${${plugin}_GIT_BRANCH} ${MAX_BRANCH_WIDTH} branch_padded)
            message(STATUS "${plugin_padded}    Branch: ${branch_padded} Revision: ${${plugin}_GIT_REVISION}")
        endif ()
    endforeach()
endfunction()

function(ListToolsGit)
    set(MAX_NAME_WIDTH 25)
    set(MAX_BRANCH_WIDTH 30)
    file(GLOB subdirs RELATIVE "${CMAKE_SOURCE_DIR}/tools" "${CMAKE_SOURCE_DIR}/tools/*")
    foreach(subdir IN LISTS subdirs)
        get_git_info(${CMAKE_SOURCE_DIR}/tools/${subdir} "${subdir}_GIT")
        pad_right(${subdir} ${MAX_NAME_WIDTH} subdir_padded)
        pad_right(${${subdir}_GIT_BRANCH} ${MAX_BRANCH_WIDTH} branch_padded)
        message(STATUS "${subdir_padded}    Branch: ${branch_padded} Revision: ${${subdir}_GIT_REVISION}")
    endforeach ()
endfunction()

function(ListAppsGit)
    set(MAX_NAME_WIDTH 25)
    set(MAX_BRANCH_WIDTH 30)
    file(GLOB subdirs RELATIVE "${CMAKE_SOURCE_DIR}/apps" "${CMAKE_SOURCE_DIR}/apps/*")
    foreach(subdir IN LISTS subdirs)
        get_git_info(${CMAKE_SOURCE_DIR}/apps/${subdir} "${subdir}_GIT")
        pad_right(${subdir} ${MAX_NAME_WIDTH} subdir_padded)
        pad_right(${${subdir}_GIT_BRANCH} ${MAX_BRANCH_WIDTH} branch_padded)
        message(STATUS "${subdir_padded}    Branch: ${branch_padded} Revision: ${${subdir}_GIT_REVISION}")
    endforeach ()
endfunction()

function(ListExternalsGit)
    set(MAX_NAME_WIDTH 25)
    set(MAX_BRANCH_WIDTH 30)
    file(GLOB subdirs RELATIVE "${CMAKE_SOURCE_DIR}/externals" "${CMAKE_SOURCE_DIR}/externals/*")
    foreach(subdir IN LISTS subdirs)
        get_git_info(${CMAKE_SOURCE_DIR}/externals/${subdir} "${subdir}_GIT")
        pad_right(${subdir} ${MAX_NAME_WIDTH} subdir_padded)
        pad_right(${${subdir}_GIT_BRANCH} ${MAX_BRANCH_WIDTH} branch_padded)
        message(STATUS "${subdir_padded}    Branch: ${branch_padded} Revision: ${${subdir}_GIT_REVISION}")
    endforeach ()
endfunction()

function(pad_right input width output_var)
    string(LENGTH "${input}" len)
    math(EXPR pad_len "${width} - ${len}")
    string(REPEAT " " ${pad_len} padding)
    set(${output_var} "${input}${padding}" PARENT_SCOPE)
endfunction()

function(ListGitRevision)
    set(MAX_NAME_WIDTH 25)
    set(MAX_BRANCH_WIDTH 30)
    message(STATUS "ø")
    message(STATUS "============ Core    ============")
    get_git_info(${CMAKE_SOURCE_DIR}/ugcore "UGCORE_GIT")
    pad_right("ugcore" ${MAX_NAME_WIDTH} ugcore_padded)
    pad_right(${UGCORE_GIT_BRANCH} ${MAX_BRANCH_WIDTH} branch_padded)
    message(STATUS "${ugcore_padded}    Branch: ${branch_padded} Revision: ${UGCORE_GIT_REVISION}")

    message(STATUS "============ Externals  =========")
    ListExternalsGit()
    message(STATUS "============ Plugins ============")
    ListPluginsGit()
    message(STATUS "============ Tools   ============")
    ListToolsGit()
    message(STATUS "============ Apps    ============")
    ListAppsGit()
    message(STATUS "---------------------------------")

endfunction()
