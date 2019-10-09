module.exports = function (grunt) {
  // Load the plugin that provides the "uglify" task.
  grunt.loadNpmTasks('grunt-contrib-uglify');
  grunt.loadNpmTasks('grunt-babel');
  grunt.loadNpmTasks('grunt-mocha-test');
  grunt.loadNpmTasks('grunt-eslint');

  grunt.initConfig({
    pkg: grunt.file.readJSON('package.json'),
    babel: {
      options: {
        sourceMap: true,
        presets: ['@babel/preset-env'],
      },
      dist: {
        files: {
          'build/jsgraphs.js': 'src/jsgraphs.js'
        },
      },
    },
    uglify: {
      options: {
        preserveComments: 'some',
      },
      build: {
        src: 'build/jsgraphs.js',
        dest: 'build/jsgraphs.min.js',
      },
    },
    mochaTest: {
      test: {
        options: {
          reporter: 'spec',
        },
        src: ['tests/**/*.js'],
      },
    },
  });

  // Default task(s).
  grunt.registerTask('default', ['babel', 'uglify']);
  grunt.registerTask('test', ['mochaTest']);
};