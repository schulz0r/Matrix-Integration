cmake_minimum_required(VERSION 3.10)
project(Matrix-Integration LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)

find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)
		 
file(GLOB SRC "*.cpp")
file(GLOB HEADERS "include/*.h")

add_executable(${PROJECT_NAME} ${SRC} ${HEADERS})

target_include_directories(${PROJECT_NAME}
	PRIVATE
	include
)

target_link_libraries(${PROJECT_NAME} PRIVATE
	Boost::boost
	Eigen3::Eigen
    )

enable_testing()

add_test(test ${PROJECT_NAME})
